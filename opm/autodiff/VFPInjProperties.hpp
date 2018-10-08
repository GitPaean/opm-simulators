/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFF_VFPINJPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPINJPROPERTIES_HPP_

#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/core/wells.h>
// TODO: this one not sure should be here.
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/autodiff/VFPHelpers.hpp>

#include <vector>
#include <map>



namespace Opm {

template <class Scalar>
class AutoDiffBlock;

class VFPInjProperties {
public:
    typedef AutoDiffBlock<double> ADB;

    /**
     * Empty constructor
     */
    VFPInjProperties();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_table A *single* VFPINJ table
     */
    explicit VFPInjProperties(const VFPInjTable* inj_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     */
    explicit VFPInjProperties(const std::map<int, std::shared_ptr<const VFPInjTable> >& inj_tables);

    /**
     * Linear interpolation of bhp as function of the input parameters.
     * @param table_id Table number to use
     * @param wells Wells structure with information about wells in qs
     * @param qs Flow quantities
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    ADB bhp(const std::vector<int>& table_id,
            const Wells& wells,
            const ADB& qs,
            const ADB& thp) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters given as ADBs
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    ADB bhp(const std::vector<int>& table_id,
            const ADB& aqua,
            const ADB& liquid,
            const ADB& vapour,
            const ADB& thp) const;

    /**
     * Linear interpolation of bhp as a function of the input parameters given as
     * Evaluation
     * Each entry corresponds typically to one well.
     * @param table_id Table number to use. A negative entry (e.g., -1)
     *                 will indicate that no table is used, and the corresponding
     *                 BHP will be calculated as a constant -1e100.
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table, for each entry in the
     * input ADB objects.
     */
    template <class EvalWell>
    EvalWell bhp(const int table_id,
                 const EvalWell& aqua,
                 const EvalWell& liquid,
                 const EvalWell& vapour,
                 const double& thp) const {

        //Get the table
        const VFPInjTable* table = detail::getTable(m_tables, table_id);
        EvalWell bhp = 0.0;

        //Find interpolation variables
        EvalWell flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());

        //Compute the BHP for each well independently
        if (table != nullptr) {
            //First, find the values to interpolate between
            //Value of FLO is negative in OPM for producers, but positive in VFP table
            auto flo_i = detail::findInterpData(flo.value(), table->getFloAxis());
            auto thp_i = detail::findInterpData( thp, table->getTHPAxis()); // assume constant

            detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i);

            bhp = bhp_val.dflo * flo;
            bhp.setValue(bhp_val.value); // thp is assumed constant i.e.
        }
        else {
            bhp.setValue(-1e100); //Signal that this value has not been calculated properly, due to "missing" table
        }
        return bhp;
    }

    /**
     * Linear interpolation of bhp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param thp Tubing head pressure
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double bhp(int table_id,
               const double& aqua,
               const double& liquid,
               const double& vapour,
               const double& thp) const;


    /**
     * Linear interpolation of thp as a function of the input parameters
     * @param table_id Table number to use
     * @param aqua Water phase
     * @param liquid Oil phase
     * @param vapour Gas phase
     * @param bhp Bottom hole pressure
     *
     * @return The tubing hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double thp(int table_id,
               const double& aqua,
               const double& liquid,
               const double& vapour,
               const double& bhp) const;

    /**
     * Returns the table associated with the ID, or throws an exception if
     * the table does not exist
     */
    const VFPInjTable* getTable(const int table_id) const;

    /**
     * Returns true if no vfp tables are in the current map
     */
    bool empty() const {
        return m_tables.empty();
    }


    /**
     * Return the bhp value based on the input flo rate
     */
    double bhp(const int table_id,
               const double& flo,
               const double& thp) const {

        //Get the table
        const VFPInjTable* table = detail::getTable(m_tables, table_id);

        //First, find the values to interpolate between
        //Value of FLO is negative in OPM for producers, but positive in VFP table
        const auto flo_i = detail::findInterpData(flo, table->getFloAxis());
        const auto thp_i = detail::findInterpData(thp, table->getTHPAxis()); // assume constant

        const detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i);

        return bhp_val.value;
    }

    // TODO: might need to be moved to cpp file
    double calculateBhpWithTHPTarget(const std::vector<double>& rates1,
                                     const std::vector<double>& rates2,
                                     // bhp2 is the bhp limit, the bhp1 is the middle af bhp2 and cell pressure
                                     const double bhp1,
                                     const double bhp2,
                                     const double dp,
                                     const int table_id,
                                     const double thp) const
    {
        assert(table_id > 0);

        const VFPInjTable* table = detail::getTable(m_tables, table_id);

        const int Water = BlackoilPhases::Aqua;
        const int Oil = BlackoilPhases::Liquid;
        const int Gas = BlackoilPhases::Vapour;

        // FLO is the rate
        const double aqua1 = rates1[Water];
        const double liquid1 = rates1[Oil];
        const double vapour1 = rates1[Gas];
        const double flo_rate1 = detail::getFlo(aqua1, liquid1, vapour1, table->getFloType());

        const double aqua2 = rates2[Water];
        const double liquid2 = rates2[Oil];
        const double vapour2 = rates2[Gas];
        const double flo_rate2 = detail::getFlo(aqua2, liquid2, vapour2, table->getFloType());

        /* std::cout << "flo_rate1 is " << flo_rate1 << std::endl;
        std::cout << " bhp1 is " << bhp1 << std::endl;

        std::cout << "flo_rate2 is " << flo_rate2 << std::endl;
        std::cout << " bhp2 is " << bhp2 << std::endl; */

        const int sample_number = 1000;
        std::vector<double> bhp_samples(sample_number);
        std::vector<double> rate_samples(sample_number, 0.);

        // rate sampling interval
        const double rate_iterval = flo_rate2 / (sample_number - 1);
        for (int i = 0; i < sample_number; ++i) {
            rate_samples[i] =  i * rate_iterval;
        }

        // based on all the rate samples, let us calculate the bhp_samples
        for (int i = 0; i < sample_number; ++i) {
            bhp_samples[i] = bhp(table_id, rate_samples[i], thp) - dp;
        }

        /* std::cout << " the rate and bhp samples " << std::endl;
        for (int i = 0; i < sample_number; ++i) {
            std::cout << rate_samples[i] << " " << bhp_samples[i] << std::endl;
        } */

        double return_bhp = 1.e-100;
        // interface needs to be resigned
        const bool found = detail::findIntersectionForBhp(rate_samples, bhp_samples, flo_rate1, flo_rate2, bhp1, bhp2, return_bhp);

        // TODO: we need to do some sanity check of the value of the return_bhp. It can be negative. Or even below the reservoir pressure
        // is not allowed, since it means it can not inject under this THP constraint. But the value of dp is not very accurate here, it can
        // result in the situation wrongly.

        // std::cout << " return_bhp " << return_bhp << std::endl;

        if (!found) {
            std::cout << " COULD NOT find an Intersection point, the well might need to be closed " << std::endl;
            std::cin.ignore();
        }

        if (return_bhp < 1.0e5) {
            std::cout << " a non-realistic bhp is obtained, we need to check it " << std::endl;
            std::cin.ignore();
        }

        return return_bhp;
    }


private:
    // Map which connects the table number with the table itself
    std::map<int, const VFPInjTable*> m_tables;
};



} //namespace



#endif /* OPM_AUTODIFF_VFPINJPROPERTIES_HPP_ */
