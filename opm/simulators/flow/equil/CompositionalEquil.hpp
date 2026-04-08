// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/**
 * \file
 *
 * \brief Compositional equilibration: replaces Black-Oil PVT density
 *        calculations with EoS flash-based density for the hydrostatic
 *        pressure ODE dP/dz = rho(z,P) * g.
 *
 * Design philosophy:
 *   - Reuse the existing RK4IVP solver, PressureFunction, and cellLoop framework
 *   - Only replace the PhasePressODE density classes with flash-based versions
 *   - Provide a CompositionalPressureTable that handles composition-vs-depth
 *     (ZMFVD) and uses a single combined ODE that tracks all phase pressures
 *     together
 */
#ifndef OPM_COMPOSITIONAL_EQUIL_HPP
#define OPM_COMPOSITIONAL_EQUIL_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/Constants.hpp>

#include <opm/simulators/flow/equil/InitStateEquil.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Opm {
namespace EQUIL {
namespace Details {

// ============================================================================
// Composition-vs-Depth interpolator
// ============================================================================

/// Interpolates component mole fractions as a function of depth.
/// Built from ZMFVD table data (one table per EQUIL region).
/// If no ZMFVD is provided, returns a constant composition.
template <class Scalar>
class CompositionVsDepth
{
public:
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    /// Construct from constant composition (no ZMFVD).
    /// \param[in] zi  Overall mole fractions for each component.
    explicit CompositionVsDepth(const std::vector<Scalar>& zi)
        : numComp_(zi.size())
        , constant_(true)
        , zi_(zi)
    {}

    /// Construct from ZMFVD table data.
    /// \param[in] depths        Depth nodes from ZMFVD.
    /// \param[in] compositions  compositions[comp][row] = z_comp at that depth.
    /// \param[in] numComp       Number of components.
    CompositionVsDepth(const std::vector<Scalar>& depths,
                       const std::vector<std::vector<Scalar>>& compositions,
                       const int numComp)
        : numComp_(numComp)
        , constant_(false)
    {
        zTables_.resize(numComp);
        for (int c = 0; c < numComp; ++c) {
            zTables_[c].setXYArrays(depths.size(), depths, compositions[c]);
        }
    }

    /// Evaluate overall mole fraction of component c at given depth.
    Scalar operator()(const int c, const Scalar depth) const
    {
        if (constant_) {
            return zi_[c];
        }
        return zTables_[c].eval(depth, /*extrapolate=*/true);
    }

    /// Return all component mole fractions at the given depth.
std::vector<Scalar> eval(const Scalar depth) const
{
    std::vector<Scalar> z(numComp_);
    Scalar sum = 0.0;

    for (int c = 0; c < numComp_; ++c) {
        const Scalar clampedDepth = std::clamp(depth,
                                            zTables_[c].xMin(),
                                            zTables_[c].xMax());
        z[c] = std::clamp((*this)(c, clampedDepth), Scalar{1.e-10}, Scalar{1.0});
        sum += z[c];
    }
    // Normalize to ensure sum = 1
    if (sum > 0.0) {
        for (int c = 0; c < numComp_; ++c) z[c] /= sum;
    }
    return z;
}

    int numComponents() const { return numComp_; }

private:
    int numComp_;
    bool constant_;
    std::vector<Scalar> zi_;                  // Used when constant_
    std::vector<TabulatedFunction> zTables_;   // Used when !constant_
};

// ============================================================================
// Flash-based density calculator
// ============================================================================

/// Computes the density of a given phase at (depth, pressure) using
/// an EoS flash calculation.  This replaces PhasePressODE::Oil/Gas/Water
/// for compositional models.
///
/// The key idea: given (T, P, z_i) at a depth, perform a PT flash to
/// determine phase split, then compute the density of the requested phase
/// from the EoS molar volume.
///
/// For the hydrostatic ODE, we need:
///   operator()(depth, press) = density(depth, press) * g
///
/// Template parameters:
///   FluidSystem  - compositional fluid system (e.g., ThreeComponentFluidSystem or GenericOil/Gas)
///   phaseIdx     - which phase pressure ODE we are solving (oil or gas)
template <class FluidSystem, unsigned phaseIdx>
class CompositionalPhasePressODE
{
    using Scalar = typename FluidSystem::Scalar;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using FlashSolver = PTFlash<Scalar, FluidSystem>;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
    using EOSType = CompositionalConfig::EOSType;
    using ParamCache = typename FluidSystem::template ParameterCache<Scalar>;

    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    /// Constructor.
    ///
    /// \param[in] tempVdTable      Temperature-vs-depth table.
    /// \param[in] compositionVsDepth  Composition interpolator (from ZMFVD or constant).
    /// \param[in] eosType          Equation of state type (PR, SRK, etc.).
    /// \param[in] normGrav         Norm of gravity vector (typically 9.80665).
    CompositionalPhasePressODE(const TabulatedFunction& tempVdTable,
                               const CompositionVsDepth<Scalar>& compositionVsDepth,
                               const EOSType eosType,
                               const Scalar normGrav)
        : tempVdTable_(tempVdTable)
        , zVsDepth_(compositionVsDepth)
        , eosType_(eosType)
        , g_(normGrav)
    {}

    /// RHS of the hydrostatic ODE: dP/dz = rho(z,P) * g.
    /// \param[in] depth  Current depth.
    /// \param[in] press  Current phase pressure estimate.
    /// \return  density * gravity
    Scalar operator()(const Scalar depth,
                      const Scalar press) const
    {
        return this->density(depth, press) * g_;
    }

private:
    const TabulatedFunction& tempVdTable_;
    const CompositionVsDepth<Scalar>& zVsDepth_;
    EOSType eosType_;
    Scalar g_;

    /// Compute the phase density at (depth, press) using an EoS flash.
    ///
    /// Steps:
    ///   1. Look up T(depth) from TEMPVD table.
    ///   2. Look up z_i(depth) from ZMFVD table (or constant composition).
    ///   3. Set up CompositionalFluidState with (T, P, z_i).
    ///   4. Perform PT flash to get phase compositions (x_i, y_i) and
    ///      phase fractions (L, V).
    ///   5. Compute phase density = averageMolarMass / molarVolume
    ///      from the EoS parameter cache.
    ///
    /// If the flash says the phase doesn't exist (e.g., single-phase region),
    /// we still return a density — either the existing phase's density
    /// or an extrapolated value.  This is analogous to how the Black-Oil
    /// code handles undersaturated oil/gas.
    Scalar density(const Scalar depth, const Scalar press) const
    {
        const Scalar temp = tempVdTable_.eval(depth, /*extrapolate=*/true);
        auto zi = zVsDepth_.eval(depth);

        // Set up a CompositionalFluidState
        FluidState fluidState;
        fluidState.setTemperature(temp);

        // Set pressures — for equilibration, we assume capillary pressure
        // between oil and gas is zero (or small) at this stage.
        // The pressure argument is for the phase whose ODE we're solving.
        fluidState.setPressure(oilPhaseIdx, press);
        fluidState.setPressure(gasPhaseIdx, press);
        if constexpr (FluidSystem::numPhases > 2) {
            // If water phase exists, set it to the same pressure for now
            fluidState.setPressure(FluidSystem::waterPhaseIdx, press);
        }

        // Set overall compositions as initial guess for flash
        for (int c = 0; c < numComponents; ++c) {
            fluidState.setMoleFraction(c, std::max(zi[c], Scalar{1.e-10}));
        }

        std::cout << " at depth " << depth << " pressure " << press << " temp " << temp << " zi " << zi[0] << " " << zi[1] << " " << zi[2] << std::endl;

        // Initialize K-values with Wilson correlation as initial guess:
        //   K_i = (Pc_i / P) * exp(5.373 * (1 + omega_i) * (1 - Tc_i/T))
        for (int c = 0; c < numComponents; ++c) {
            const Scalar Tc = FluidSystem::criticalTemperature(c);
            const Scalar Pc = FluidSystem::criticalPressure(c);
             const Scalar omega = FluidSystem::acentricFactor(c);
            Scalar Ki = (Pc / press) * std::exp(5.373 * (1.0 + omega) * (1.0 - Tc / temp));
            Ki = std::max(Ki, Scalar{1e-10});
            Ki = std::min(Ki, Scalar{1e10});
            fluidState.setKvalue(c, Ki);
        }

        // Initial guess for liquid fraction using Rachford-Rice perspective
        fluidState.setLvalue(Scalar{0.5});

        // Perform PT flash (scalar version — no AD derivatives)
        const bool isSinglePhase = FlashSolver::flash_solve_scalar_(
            fluidState, "ssi", /*tolerance=*/1e-7, eosType_, /*verbosity=*/0
        );

        // Compute the molar volume and density of the requested phase
        // using the EoS parameter cache.
        ParamCache paramCache(eosType_);

        Scalar density {0.};

        if (isSinglePhase) {
            // Single phase: use whichever phase exists.
            // If L ~ 1 => only liquid (oil), if L ~ 0 => only vapor (gas).
            // TODO: using the approach in CompWell_impl.hpp to do the single phase labeling
            const Scalar L = fluidState.L();
            const unsigned existingPhase = (L > 0.5) ? oilPhaseIdx : gasPhaseIdx;

            // Even though our ODE is for `phaseIdx`, in a single-phase
            // region, the single existing phase fills the entire pore space.
            // We return its density regardless of which ODE we're solving,
            // since the pressure of a non-existing phase is extrapolated.
            paramCache.updatePhase(fluidState, existingPhase);
            const Scalar Vm = paramCache.molarVolume(existingPhase);
            const Scalar avgMW = fluidState.averageMolarMass(existingPhase);
            density = avgMW / Vm;
        }
        else {
            // TODO: we need to think about how the density should be calculated
            // Two-phase: compute density of the specific phase
            paramCache.updatePhase(fluidState, phaseIdx);
            const Scalar Vm = paramCache.molarVolume(phaseIdx);
            const Scalar avgMW = fluidState.averageMolarMass(phaseIdx);
            density = avgMW / Vm;
        }
        std::cout << " is single phase ? " << isSinglePhase << " density of phase " << phaseIdx << " is " << density << std::endl;
        return density;
    }
};

// ============================================================================
// Compositional PressureTable
// ============================================================================

/// Manages phase pressure-vs-depth tables for compositional models.
///
/// Analogous to Details::PressureTable<FluidSystem, Region> for Black-Oil,
/// but uses CompositionalPhasePressODE instead of PhasePressODE::Oil/Gas/Water.
///
/// The equilibration strategy is simpler for compositional models:
///   - Gas and oil pressures are both computed from the datum,
///     with P_gas - P_oil = 0 at the datum (or specified Pcgo).
///   - Water phase (if present) uses the Black-Oil water PVT as before,
///     since water is typically not part of the compositional system.
///
/// Template parameters:
///   FluidSystem - compositional fluid system
///   Region      - equilibration region type (EquilReg<Scalar>)
template <class FluidSystem, class Region>
class CompositionalPressureTable
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using VSpan = std::array<Scalar, 2>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using EOSType = CompositionalConfig::EOSType;

    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

private:
    // Phase pressure ODE types
    using OilPressODE = CompositionalPhasePressODE<FluidSystem, oilPhaseIdx>;
    using GasPressODE = CompositionalPhasePressODE<FluidSystem, gasPhaseIdx>;

    // Reuse PressureFunction from InitStateEquil (it's a generic wrapper
    // around any ODE + RK4IVP).
    template <class ODE>
    class PressureFunction
    {
    public:
        struct InitCond {
            Scalar depth;
            Scalar pressure;
        };

        PressureFunction(const ODE&      ode,
                         const InitCond& ic,
                         const int       nsample,
                         const VSpan&    span)
            : initial_(ic)
        {
            value_[Direction::Up] = std::make_unique<Distribution>(
                ode, VSpan{{ ic.depth, span[0] }}, ic.pressure, nsample);
            value_[Direction::Down] = std::make_unique<Distribution>(
                ode, VSpan{{ ic.depth, span[1] }}, ic.pressure, nsample);
        }

        PressureFunction(const PressureFunction& rhs)
            : initial_(rhs.initial_)
        {
            value_[Direction::Up]   = std::make_unique<Distribution>(*rhs.value_[Direction::Up]);
            value_[Direction::Down] = std::make_unique<Distribution>(*rhs.value_[Direction::Down]);
        }

        PressureFunction(PressureFunction&& rhs) = default;

        PressureFunction& operator=(const PressureFunction& rhs)
        {
            initial_ = rhs.initial_;
            value_[Direction::Up]   = std::make_unique<Distribution>(*rhs.value_[Direction::Up]);
            value_[Direction::Down] = std::make_unique<Distribution>(*rhs.value_[Direction::Down]);
            return *this;
        }

        PressureFunction& operator=(PressureFunction&&) = default;

        Scalar value(const Scalar depth) const
        {
            if (depth < initial_.depth)
                return (*value_[Direction::Up])(depth);
            else if (depth > initial_.depth)
                return (*value_[Direction::Down])(depth);
            else
                return initial_.pressure;
        }

    private:
        enum Direction : std::size_t { Up, Down, NumDir };
        using Distribution = RK4IVP<Scalar, ODE>;
        using DistrPtr = std::unique_ptr<Distribution>;

        InitCond initial_;
        std::array<DistrPtr, Direction::NumDir> value_;
    };

    using OPress = PressureFunction<OilPressODE>;
    using GPress = PressureFunction<GasPressODE>;

public:
    /// Constructor.
    ///
    /// \param[in] gravity       Gravitational acceleration (m/s^2).
    /// \param[in] samplePoints  Number of RK4 sample points per table.
    /// \param[in] compositionVsDepth  Composition interpolator.
    /// \param[in] eosType       Equation of state type.
    CompositionalPressureTable(const Scalar gravity,
                               const int    samplePoints,
                               const CompositionVsDepth<Scalar>& compositionVsDepth,
                               const EOSType eosType)
        : gravity_(gravity)
        , nsample_(samplePoints)
        , zVsDepth_(compositionVsDepth)
        , eosType_(eosType)
    {}

    /// Compute equilibration pressure tables for all phases.
    ///
    /// For compositional models the strategy is:
    ///   1. The datum depth gives a reference pressure for the oil phase.
    ///   2. Solve dP_oil/dz = rho_oil(z, P_oil) * g from the datum depth.
    ///   3. At the GOC, P_gas = P_oil + Pcgo.
    ///   4. Solve dP_gas/dz = rho_gas(z, P_gas) * g from the GOC.
    ///   5. Water (if present) handled as in Black-Oil.
    ///
    /// For simplicity, in the initial implementation we assume:
    ///   - The datum is in the oil zone.
    ///   - If no distinct GOC, gas and oil pressure are computed from
    ///     the same datum with the same reference pressure (Pcgo=0).
    void equilibrate(const Region& reg,
                     const VSpan&  span)
    {
        const auto& tempVdTable = reg.tempVdTable();

        // TODO: this EQUIL setup, the datum depth is in the gas region, we should do the gas region first
        // same with the blackoil version, we should do the different strategy
        // then the gas-oil contact point, the gas pressure is the same as oil pressure,
        // then based on this pressure, we can calculate the oil pressure

        // --- gas pressure ---
        {
            const auto gasODE = GasPressODE{
                tempVdTable, zVsDepth_, eosType_, gravity_
            };
            const auto ic = typename GPress::InitCond{
                reg.datum(), reg.pressure()
            };
            gas_ = std::make_unique<GPress>(gasODE, ic, nsample_, span);
        }

        // --- oil pressure ---
        {
            const auto oilODE = OilPressODE{
                tempVdTable, zVsDepth_, eosType_, gravity_
            };

            // Gas initial condition from GOC:
            //   P_oil(GOC) = P_gas(GOC) - Pcgo
            const Scalar zgoc = reg.zgoc();
            const Scalar pGasAtGoc = gas_->value(zgoc);
            const Scalar pcgoGoc = reg.pcgoGoc();

            const auto ic = typename OPress::InitCond{
                zgoc, pGasAtGoc - pcgoGoc
            };
            oil_ = std::make_unique<OPress>(oilODE, ic, nsample_, span);
        }

        // --- Water pressure ---
        // For compositional models, water is typically treated as a separate
        // immiscible phase with Black-Oil water PVT.  If we have water,
        // we compute it using the standard approach:
        //   P_wat(WOC) = P_oil(WOC) - Pcow
        // This part is left to the existing Black-Oil water pressure machinery.
        // The caller should set up the water pressure table separately
        // if water is active.
    }

    /// Evaluate oil phase pressure at specified depth.
    Scalar oil(const Scalar depth) const
    {
        if (!oil_) {
            throw std::logic_error("CompositionalPressureTable: oil pressure not initialized");
        }
        return oil_->value(depth);
    }

    /// Evaluate gas phase pressure at specified depth.
    Scalar gas(const Scalar depth) const
    {
        if (!gas_) {
            throw std::logic_error("CompositionalPressureTable: gas pressure not initialized");
        }
        return gas_->value(depth);
    }

    bool oilActive() const { return oil_ != nullptr; }
    bool gasActive() const { return gas_ != nullptr; }

    /// Composition at a given depth (for post-processing).
    std::vector<Scalar> composition(const Scalar depth) const
    {
        return zVsDepth_.eval(depth);
    }

    const CompositionVsDepth<Scalar>& compositionVsDepth() const
    {
        return zVsDepth_;
    }

private:
    Scalar gravity_;
    int nsample_;
    CompositionVsDepth<Scalar> zVsDepth_;
    EOSType eosType_;

    std::unique_ptr<OPress> oil_{};
    std::unique_ptr<GPress> gas_{};
};

// ============================================================================
// Compositional phase saturation calculator
// ============================================================================

/// Determines phase saturations for compositional models.
///
/// In compositional equilibration, saturations come from the flash
/// calculation rather than from inverting capillary pressure curves.
/// Given the phase pressures from the CompositionalPressureTable and
/// the composition from ZMFVD, a flash at each cell centre gives
/// the phase split (L, V), from which saturations are derived.
// TODO: We should not do this because we should use the capillary pressure to do the saturation calculaiton
template <class FluidSystem>
class CompositionalPhaseSaturations
{
    using Scalar = typename FluidSystem::Scalar;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
    using FlashSolver = PTFlash<Scalar, FluidSystem>;
    using EOSType = CompositionalConfig::EOSType;
    using ParamCache = typename FluidSystem::template ParameterCache<Scalar>;

    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    struct Result {
        Scalar So;  // oil saturation
        Scalar Sg;  // gas saturation
        Scalar Sw;  // water saturation (if applicable)
        Scalar rho_oil;
        Scalar rho_gas;
        std::vector<Scalar> x;   // oil phase mole fractions [numComponents]
        std::vector<Scalar> y;   // gas phase mole fractions [numComponents]
        Scalar L;   // liquid (oil) molar fraction
    };

    explicit CompositionalPhaseSaturations(const EOSType eosType)
        : eosType_(eosType)
    {}

    /// Compute phase saturations at (depth, P_oil, P_gas, T, z_i).
    ///
    /// This performs a flash calculation and converts the molar phase
    /// split into volume-based saturations:
    ///
    ///   S_oil = (L / rho_oil_molar) / (L/rho_oil_molar + V/rho_gas_molar)
    ///   S_gas = 1 - S_oil   (for two-phase oil+gas)
    ///
    /// where rho_*_molar = 1 / Vm_* is the molar density from the EoS.
    Result compute(const Scalar pOil,
                   const Scalar pGas,
                   const Scalar temp,
                   const std::vector<Scalar>& zi) const
    {
        Result result{};
        result.x.resize(numComponents);
        result.y.resize(numComponents);

        FluidState fluidState;
        fluidState.setTemperature(temp);
        fluidState.setPressure(oilPhaseIdx, pOil);
        fluidState.setPressure(gasPhaseIdx, pGas);

        // Set overall mole fractions
        for (int c = 0; c < numComponents; ++c) {
            fluidState.setMoleFraction(c, std::max(zi[c], 1.e-10));
        }

        // Wilson K-value initial guess
        for (int c = 0; c < numComponents; ++c) {
            const Scalar Tc = FluidSystem::criticalTemperature(c);
            const Scalar Pc = FluidSystem::criticalPressure(c);
            const Scalar omega = FluidSystem::acentricFactor(c);
            Scalar Ki = (Pc / pOil) * std::exp(5.373 * (1.0 + omega) * (1.0 - Tc / temp));
            Ki = std::clamp(Ki, Scalar{1e-10}, Scalar{1e10});
            fluidState.setKvalue(c, Ki);
        }
        fluidState.setLvalue(Scalar{0.5});

        const bool isSinglePhase = FlashSolver::flash_solve_scalar_(
            fluidState, "ssi", 1e-7, eosType_, 0
        );

        ParamCache paramCache(eosType_);

        if (isSinglePhase) {
            // TODO: we should the slightly better way in CompWell_impl.hpp to mark the phase
            // and calculate the saturations
            const Scalar L = fluidState.L();
            if (L > 0.5) {
                // All liquid (oil)
                result.So = 1.0;
                result.Sg = 0.0;
                result.L  = 1.0;
            } else {
                // All vapor (gas)
                result.So = 0.0;
                result.Sg = 1.0;
                result.L  = 0.0;
            }
        }
        else {
            // Two-phase: convert molar split to volumetric saturations.
            //
            //   L = moles of liquid / total moles
            //   V = 1 - L
            //
            //   Volume_oil = L * Vm_oil
            //   Volume_gas = V * Vm_gas
            //
            //   S_oil = Volume_oil / (Volume_oil + Volume_gas)
            //   S_gas = 1 - S_oil
            const Scalar L = fluidState.L();
            const Scalar V = 1.0 - L;

            paramCache.updatePhase(fluidState, oilPhaseIdx);
            const Scalar VmOil = paramCache.molarVolume(oilPhaseIdx);

            paramCache.updatePhase(fluidState, gasPhaseIdx);
            const Scalar VmGas = paramCache.molarVolume(gasPhaseIdx);

            const Scalar volOil = L * VmOil;
            const Scalar volGas = V * VmGas;
            const Scalar totalVol = volOil + volGas;

            result.So = volOil / totalVol;
            result.Sg = volGas / totalVol;
            result.L  = L;
        }

        // Store phase compositions
        for (int c = 0; c < numComponents; ++c) {
            result.x[c] = fluidState.moleFraction(oilPhaseIdx, c);
            result.y[c] = fluidState.moleFraction(gasPhaseIdx, c);
        }

        // Compute densities
        if (result.So > 0.0) {
            paramCache.updatePhase(fluidState, oilPhaseIdx);
            result.rho_oil = fluidState.averageMolarMass(oilPhaseIdx) /
                             paramCache.molarVolume(oilPhaseIdx);
        }
        if (result.Sg > 0.0) {
            paramCache.updatePhase(fluidState, gasPhaseIdx);
            result.rho_gas = fluidState.averageMolarMass(gasPhaseIdx) /
                             paramCache.molarVolume(gasPhaseIdx);
        }

        result.Sw = 0.0; // Water handled separately
        return result;
    }

private:
    EOSType eosType_;
};

// ============================================================================
// Integration helper: compositional equilibrateCellCentres
// ============================================================================

/// This function template shows how the compositional equilibration
/// integrates with the existing cellLoop framework.
///
/// It is meant to be called from InitialStateComputer::calcPressSatRsRv()
/// (or a compositional equivalent) instead of the Black-Oil
/// equilibrateCellCentres.
///
/// Usage sketch (inside calcPressSatRsRv or similar):
///
///   auto compPtable = CompositionalPressureTable<FluidSystem, EquilReg<Scalar>>{
///       grav, num_pressure_points_, compositionVsDepth, eosType
///   };
///   compPtable.equilibrate(eqreg, vspan);
///
///   auto compPhaseSat = CompositionalPhaseSaturations<FluidSystem>{ eosType };
///
///   compositionalEquilibrateCellCentres(
///       cells, eqreg, compPtable, compPhaseSat, compositionVsDepth,
///       temperature_, cellCenterDepth_,
///       pp_, sat_, cellCompositions_
///   );
///
template <class FluidSystem,
          class CellRange,
          class Region,
          class Scalar = typename FluidSystem::Scalar>
void compositionalEquilibrateCellCentres(
    const CellRange& cells,
    [[maybe_unused]] const Region& eqreg,
    const CompositionalPressureTable<FluidSystem, Region>& ptable,
    const CompositionalPhaseSaturations<FluidSystem>& phaseSatCalc,
    const CompositionVsDepth<Scalar>& zVsDepth,
    const std::vector<Scalar>& temperature,
    const std::vector<Scalar>& cellCenterDepth,
    // Output arrays:
    std::vector<std::vector<Scalar>>& pp,        // pp[phaseIdx][cell]
    std::vector<std::vector<Scalar>>& sat,       // sat[phaseIdx][cell]
    std::vector<std::vector<Scalar>>& cellCompositions)  // cellCompositions[cell * numComp + c]
{
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int numComponents = FluidSystem::numComponents;

    for (const auto& cell : cells) {
        const Scalar depth = cellCenterDepth[cell];
        const Scalar temp  = temperature[cell];

        // Look up pressures from the precomputed tables
        const Scalar pOil = ptable.oil(depth);
        const Scalar pGas = ptable.gas(depth);

        // Look up composition at this depth
        const auto zi = zVsDepth.eval(depth);

        // Flash to get saturations
        const auto result = phaseSatCalc.compute(pOil, pGas, temp, zi);

        // Store pressures
        pp[oilPhaseIdx][cell] = pOil;
        pp[gasPhaseIdx][cell] = pGas;

        // Store saturations
        sat[oilPhaseIdx][cell] = result.So;
        sat[gasPhaseIdx][cell] = result.Sg;

        // Store per-cell compositions (flattened: cell * numComp + c)
        for (int c = 0; c < numComponents; ++c) {
            cellCompositions[cell * numComponents + c] = zi[c];
        }
    }
}

} // namespace Details
} // namespace EQUIL
} // namespace Opm

#endif // OPM_COMPOSITIONAL_EQUIL_HPP
