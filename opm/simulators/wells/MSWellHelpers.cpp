/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2020 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <opm/simulators/wells/MSWellHelpers.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif // HAVE_UMFPACK

#include <cmath>
#include <dune/istl/matrixmarket.hh>

namespace {

template <typename ValueType>
ValueType haalandFormular(const ValueType& re,
                          const double diameter,
                          const double roughness)
{
    const ValueType value = -3.6 * log10(6.9 / re + std::pow(roughness / (3.7 * diameter), 10. / 9.) );

    // sqrt(1/f) should be non-positive
    assert(value >= 0.0);

    return 1. / (value * value);
}

// water in oil emulsion viscosity
// TODO: maybe it should be two different ValueTypes. When we calculate the viscosity for transitional zone
template <typename ValueType>
    ValueType WIOEmulsionViscosity(const ValueType& oil_viscosity, const ValueType& water_liquid_fraction,
                                   const double max_visco_ratio, const bool output=false)
{
    const ValueType temp_value = 1. / (1. - (0.8415 / 0.7480 * water_liquid_fraction) );
    const ValueType viscosity_ratio = pow(temp_value, 2.5);
    if (output) {
        std::cout << " WIOEmulsionViscosity temp_value ";
        temp_value.print();
        std::cout << std::endl;
        std::cout << " viscosity_ratio ";
        viscosity_ratio.print();
        std::cout << std::endl;
        std::cout << " oil_viscosity ";
        oil_viscosity.print();
        std::cout << std::endl;
    }

    ValueType result;
    if (viscosity_ratio <= max_visco_ratio) {
        result = oil_viscosity * viscosity_ratio;
    } else {
        result = oil_viscosity * max_visco_ratio;
    }
    if (output) {
        std::cout << "WIOEmulsionViscosity result ";
        result.print();
        std::cout << std::endl;
    }
    return result;
}

// oil in water emulsion viscosity
template <typename ValueType>
    ValueType OIWEmulsionViscosity(const ValueType& water_viscosity, const ValueType& water_liquid_fraction,
                                   const double max_visco_ratio, const bool output=false)
{
    const ValueType temp_value = 1. / (1. - (0.6019 / 0.6410) * (1. - water_liquid_fraction) );
    const ValueType viscosity_ratio = pow(temp_value, 2.5);
    if (output) {
        std::cout << " OIWEmulsionViscosity temp_value ";
        temp_value.print();
        std::cout << std::endl;
        std::cout << " viscosity_ratio ";
        viscosity_ratio.print();
        std::cout << std::endl;
        std::cout << " water_viscosity ";
        water_viscosity.print();
        std::cout << std::endl;
    }
    ValueType result;
    if (viscosity_ratio <= max_visco_ratio) {
        result = water_viscosity * viscosity_ratio;
    } else {
        result = water_viscosity * max_visco_ratio;
    }
    if (output) {
        std::cout << "OIWEmulsionViscosity result ";
        result.print();
        std::cout << std::endl;
    }
    return result;
}

}

namespace Opm {

namespace mswellhelpers
{

    /// Applies umfpack and checks for singularity
template <typename MatrixType, typename VectorType>
    VectorType
applyUMFPack(const MatrixType& D, std::shared_ptr<Dune::UMFPack<MatrixType>>& linsolver, VectorType x,
             const std::string& well_name)
{
#if HAVE_UMFPACK
    std::cout << " applyUMFPack with well " << well_name << std::endl;
    const int number = rand()%100000;
    if (!linsolver)
    {
        linsolver = std::make_shared<Dune::UMFPack<MatrixType>>(D, 0);

        const double* UMF_Decomposition_Info = linsolver->getDecompositionInfo();
        if ( (1./abs(UMF_Decomposition_Info[UMFPACK_RCOND]) > 1.e50) || (abs(UMF_Decomposition_Info[UMFPACK_RCOND]) > 1.e50)
             || !std::isfinite(UMF_Decomposition_Info[UMFPACK_RCOND])) {
            const std::string matrix_filename = "debug_output/duneD_" + std::to_string(number);
            const std::string rhs_filename = "debug_output/resWell_" + std::to_string(number);
            std::cout << "condition number is " << UMF_Decomposition_Info[UMFPACK_RCOND] << std::endl;
            std::cout << " outputting the D and rhs for well " << well_name << " to " << matrix_filename << " and " << rhs_filename << " respectively " << std::endl;
            Dune::storeMatrixMarket(D, matrix_filename);
            Dune::storeMatrixMarket(x, rhs_filename);
        }
    }

    // The copy of x seems mandatory for calling UMFPack!
    VectorType y(x.size());
    y = 0.;

    // Object storing some statistics about the solving process
    Dune::InverseOperatorResult res;

    // Solve
    linsolver->apply(y, x, res);

    // Checking if there is any inf or nan in y
    // it will be the solution before we find a way to catch the singularity of the matrix
    for (size_t i_block = 0; i_block < y.size(); ++i_block) {
        for (size_t i_elem = 0; i_elem < y[i_block].size(); ++i_elem) {
            if (std::isinf(y[i_block][i_elem]) || std::isnan(y[i_block][i_elem]) ) {
                const std::string file_name = "debug_output/solution_" + std::to_string(nubmer)+"_"+well_name;
                Dune::storeMatrixMarket(y, file_name);
                const std::string msg = "nan or inf value found after UMFPack solve due to singular matrix of well " + well_name;
                OpmLog::debug(msg);
                OPM_THROW_NOLOG(NumericalIssue, msg);
            }
        }
    }
    std::cout << " leaving applyUMFPack with well " << well_name << std::endl;
    return y;
#else
    // this is not thread safe
    OPM_THROW(std::runtime_error, "Cannot use applyUMFPack() without UMFPACK. "
              "Reconfigure opm-simulators with SuiteSparse/UMFPACK support and recompile.");
#endif // HAVE_UMFPACK
}

template <typename MatrixType, typename VectorType>
Dune::Matrix<typename MatrixType::block_type>
invertWithUMFPack(const MatrixType& D, std::shared_ptr<Dune::UMFPack<MatrixType> >& linsolver,
                  const std::string& well_name)
{
#if HAVE_UMFPACK
    const int sz = D.M();
    const int bsz = D[0][0].M();
    VectorType e(sz);
    e = 0.0;

    // Make a full block matrix.
    Dune::Matrix<typename MatrixType::block_type> inv(sz, sz);

    // Create inverse by passing basis vectors to the solver.
    for (int ii = 0; ii < sz; ++ii) {
        for (int jj = 0; jj < bsz; ++jj) {
            e[ii][jj] = 1.0;
            auto col = applyUMFPack(D, linsolver, e, well_name);
            for (int cc = 0; cc < sz; ++cc) {
                for (int dd = 0; dd < bsz; ++dd) {
                    inv[cc][ii][dd][jj] = col[cc][dd];
                }
            }
            e[ii][jj] = 0.0;
        }
    }

    return inv;
#else
    // this is not thread safe
    OPM_THROW(std::runtime_error, "Cannot use invertWithUMFPack() without UMFPACK. "
              "Reconfigure opm-simulators with SuiteSparse/UMFPACK support and recompile.");
#endif // HAVE_UMFPACK
}

template <typename MatrixType, typename VectorType>
VectorType
invDX(const MatrixType& D, VectorType x, DeferredLogger& deferred_logger)
{
    // the function will change the value of x, so we should not use reference of x here.

    // TODO: store some of the following information to avoid to call it again and again for
    // efficiency improvement.
    // Bassically, only the solve / apply step is different.

    VectorType y(x.size());
    y = 0.;

    Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(D);

    // Sequential incomplete LU decomposition as the preconditioner
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    Dune::SeqILU<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);
#else
    Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);
#endif
    // Dune::SeqILUn<MatrixType, VectorType, VectorType> preconditioner(D, 1, 0.92);
    // Dune::SeqGS<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);
    // Dune::SeqJac<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);

    // Preconditioned BICGSTAB solver
    Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
                                               preconditioner,
                                               1.e-8, // desired residual reduction factor
                                               250, // maximum number of iterations
                                               0); // verbosity of the solver */

    // Object storing some statistics about the solving process
    Dune::InverseOperatorResult res;

    // Solve
    linsolver.apply(y, x, res);

    if ( !res.converged ) {
        OPM_DEFLOG_THROW(NumericalIssue, "the invDX did not converge ", deferred_logger);
    }

    return y;
}

template <typename ValueType>
ValueType frictionPressureLoss(const double l, const double diameter,
                               const double area, const double roughness,
                               const ValueType& density,
                               const ValueType& w, const ValueType& mu)
{
    // Reynolds number
    const ValueType re = abs( diameter * w / (area * mu));

    /* if ( abs(w) <= 1.e-30 ) {
        // make sure it is because the mass rate is zero
        // assert(w == 0.);
        return 0.0;
    } */

    constexpr double re_value1 = 2000.;
    constexpr double re_value2 = 4000.;

    ValueType result;

    if (re < re_value1) {
        // not using the formula directly because of the division with w might cause singular Jacobian matrix
        result = 32.* mu * l * abs(w) / (area * diameter *diameter * density);
    } else {
        ValueType f;
        if (re > re_value2){
            f = haalandFormular(re, diameter, roughness);
        } else { // in between
            const ValueType f1 = 16. / re_value1;
            const ValueType f2 = haalandFormular(re_value2, diameter, roughness);
            f = (f2 - f1) / (re_value2 - re_value1) * (re - re_value1) + f1;
        }
        // \Note: a factor of 2 needs to be here based on the dimensional analysis
        result = 2. * f * l * w * w / (area * area * diameter * density);
    }
    return result;
}

template <typename ValueType>
ValueType valveContrictionPressureLoss(const ValueType& mass_rate,
                                       const ValueType& density,
                                       const double area_con, const double cv)
{
    // the formulation is adjusted a little bit for convinience
    // velocity = mass_rate / (density * area) is applied to the original formulation
    const double area = (area_con > 1.e-10 ? area_con : 1.e-10);
    return mass_rate * mass_rate / (2. * density * cv * cv * area * area);
}

template <typename ValueType>
ValueType velocityHead(const double area, const ValueType& mass_rate,
                       const ValueType& density)
{
    // \Note: a factor of 2 is added to the formulation in order to match results from the
    // reference simulator. This is inline with what is done for the friction loss.
    return (mass_rate * mass_rate / (area * area * density));
}

template <typename ValueType>
ValueType
emulsionViscosity(const ValueType& water_fraction, const ValueType& water_viscosity, const ValueType& oil_fraction,
                  const ValueType& oil_viscosity, const SICD& sicd, const bool output)
{
    const double width_transition = sicd.widthTransitionRegion();

    // it is just for now, we should be able to treat it.
    if (width_transition <= 0.) {
        OPM_THROW(std::runtime_error, "Not handling non-positive transition width now");
    }

    const double critical_value = sicd.criticalValue();
    const ValueType transition_start_value = critical_value - width_transition / 2.0;
    const ValueType transition_end_value = critical_value + width_transition / 2.0;

    const ValueType liquid_fraction = water_fraction + oil_fraction;
    // if there is no liquid, we just return zero
    if (liquid_fraction == 0.) {
        return 0.;
    }

    const ValueType water_liquid_fraction = water_fraction / liquid_fraction;
    if (output) {
        std::cout << " liquid_fraction ";
        liquid_fraction.print();
        std::cout << std::endl;
        std::cout << " water_liquid_fracton ";
        water_liquid_fraction.print();
        std::cout << std::endl;
    }

    const double max_visco_ratio = sicd.maxViscosityRatio();
    if (water_liquid_fraction <= transition_start_value) {
        return WIOEmulsionViscosity(oil_viscosity, water_liquid_fraction, max_visco_ratio, output);
    } else if (water_liquid_fraction >= transition_end_value) {
        return OIWEmulsionViscosity(water_viscosity, water_liquid_fraction, max_visco_ratio, output);
    } else { // in the transition region
        const ValueType viscosity_start_transition = WIOEmulsionViscosity(oil_viscosity, transition_start_value,
                                                                          max_visco_ratio, output);
        const ValueType viscosity_end_transition = OIWEmulsionViscosity(water_viscosity, transition_end_value,
                                                                        max_visco_ratio, output);
        const ValueType emulsion_viscosity = (viscosity_start_transition * (transition_end_value - water_liquid_fraction)
                                           + viscosity_end_transition * (water_liquid_fraction - transition_start_value) ) / width_transition;
        return emulsion_viscosity;
    }
}

template<int Dim>
using Vec = Dune::BlockVector<Dune::FieldVector<double,Dim>>;
template<int Dim>
using Mat = Dune::BCRSMatrix<Dune::FieldMatrix<double,Dim,Dim>>;

#define INSTANCE_UMF(Dim) \
    template Vec<Dim> applyUMFPack<Mat<Dim>,Vec<Dim>>(const Mat<Dim>&, \
                                                      std::shared_ptr<Dune::UMFPack<Mat<Dim>>>&, \
                                                      Vec<Dim>, const std::string&); \
    template Dune::Matrix<typename Mat<Dim>::block_type> \
    invertWithUMFPack<Mat<Dim>,Vec<Dim>>(const Mat<Dim>& D, \
                                         std::shared_ptr<Dune::UMFPack<Mat<Dim>>>&, const std::string&);

INSTANCE_UMF(2)
INSTANCE_UMF(3)
INSTANCE_UMF(4)

#define INSTANCE_IMPL(...) \
    template __VA_ARGS__ \
    frictionPressureLoss<__VA_ARGS__>(const double, \
                                      const double, \
                                      const double, \
                                      const double, \
                                      const __VA_ARGS__&, \
                                      const __VA_ARGS__&, \
                                      const __VA_ARGS__&); \
    template  __VA_ARGS__ \
    valveContrictionPressureLoss<__VA_ARGS__>(const __VA_ARGS__& mass_rate, \
                                              const __VA_ARGS__& density, \
                                              const double, const double); \
    template __VA_ARGS__ \
    velocityHead<__VA_ARGS__>(const double, const __VA_ARGS__&, const __VA_ARGS__&); \
    template __VA_ARGS__ \
    emulsionViscosity<__VA_ARGS__>(const __VA_ARGS__&,  \
                                   const __VA_ARGS__&, \
                                   const __VA_ARGS__&, \
                                   const __VA_ARGS__&, \
                                   const SICD&);

#define INSTANCE_EVAL(Dim) \
    INSTANCE_IMPL(DenseAd::Evaluation<double,Dim>)

INSTANCE_EVAL(3)
INSTANCE_EVAL(4)
INSTANCE_EVAL(5)
INSTANCE_EVAL(6)
INSTANCE_EVAL(7)
INSTANCE_EVAL(8)
INSTANCE_EVAL(9)

} // namespace mswellhelpers

} // namespace Opm
