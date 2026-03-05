/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS

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

#ifndef OPM_BLACKOIL_FLUID_STATE_UTILS_HPP
#define OPM_BLACKOIL_FLUID_STATE_UTILS_HPP

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/thermal/EnergyModuleType.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <fmt/format.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace Opm {

/// Creates a BlackOilFluidState from surface-volume fluid compositions.
///
/// Given per-component surface-volume fractions (\p fluid_composition), a
/// wellbore/cell pressure and temperature, a salt concentration, and a PVT
/// region index, this function computes the dissolution/vaporization
/// factors (Rs, Rv, Rsw, Rvw), the inverse formation volume factors (invB),
/// phase saturations and densities (and, when thermal is active, enthalpies),
/// and returns the fully-populated BlackOilFluidState.
///
/// This is a general-purpose factory function that is not tied to any
/// specific well or reservoir class — it only requires a TypeTag that
/// exposes FluidSystem, Indices, and the standard BlackOil property flags.
///
/// \tparam TypeTag  OPM TypeTag that provides FluidSystem, Indices,
///                  EnableBrine, EnergyModuleType, EnableVapwat, etc.
/// \tparam ValueType  Scalar or Evaluation (AD) type for the computation.
///
/// \param[in] fluid_composition  Surface-volume fractions, indexed by
///                               active component index.
/// \param[in] pressure           Phase pressure (same for all phases;
///                               capillary pressure is neglected).
/// \param[in] temperature        Fluid temperature (only used when the
///                               energy equation is active).
/// \param[in] saltConcentration  Brine salt concentration (only used when
///                               brine support is enabled).
/// \param[in] pvtRegionIdx       PVT region index for fluid property
///                               table lookups.
/// \param[in] name               Descriptive name used in error messages
///                               (e.g. well name, cell identifier).
///
/// \returns A BlackOilFluidState with pressures, dissolution factors,
///          invB, saturations, densities (and enthalpies) populated.
template <typename TypeTag, typename ValueType>
auto createBlackOilFluidState(const std::vector<ValueType>& fluid_composition,
                               const ValueType&              pressure,
                               const ValueType&              temperature,
                               GetPropType<TypeTag, Properties::Scalar> saltConcentration,
                               int         pvtRegionIdx,
                               std::string_view name = "")
{
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices     = GetPropType<TypeTag, Properties::Indices>;

    static constexpr EnergyModules energyModuleType =
        getPropValue<TypeTag, Properties::EnergyModuleType>();
    static constexpr bool has_energy =
        energyModuleType == EnergyModules::FullyImplicitThermal;
    static constexpr bool has_brine =
        getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool has_watVapor =
        getPropValue<TypeTag, Properties::EnableVapwat>();
    static constexpr bool has_saltPrecip =
        getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr bool has_disgas_in_water =
        getPropValue<TypeTag, Properties::EnableDisgasInWater>();

    using FluidState = BlackOilFluidState<ValueType,
                                          FluidSystem,
                                          energyModuleType != EnergyModules::NoTemperature,
                                          energyModuleType == EnergyModules::FullyImplicitThermal,
                                          Indices::compositionSwitchIdx >= 0,
                                          has_watVapor,
                                          has_brine,
                                          has_saltPrecip,
                                          has_disgas_in_water,
                                          Indices::numPhases>;

    FluidState fluid_state;

    if constexpr (has_energy) {
        fluid_state.setTemperature(temperature);
    }
    if constexpr (has_brine) {
        fluid_state.setSaltConcentration(saltConcentration);
    }

    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }
        // we assume there is no capillary pressure in the wellbore
        fluid_state.setPressure(phaseIdx, pressure);
    }
    fluid_state.setPvtRegionIndex(pvtRegionIdx);

    const bool both_oil_gas =
        FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
        FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const bool both_water_gas =
        FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
        FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);

    const ValueType zero_value {0.};
    // let us handle the dissolution first
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned activeCompIdx =
            FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));
        switch (phaseIdx) {
            case FluidSystem::oilPhaseIdx: {
                // TODO: we should use something else
                if constexpr (Indices::compositionSwitchIdx >= 0) {
                    if (both_oil_gas) {
                        const ValueType saturated_rs =
                            FluidSystem::saturatedDissolutionFactor(fluid_state, phaseIdx,
                                                                    fluid_state.pvtRegionIndex());
                        const unsigned gasCompIdx = FluidSystem::canonicalToActiveCompIdx(
                            FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                        const ValueType max_possible_rs =
                            fluid_composition[gasCompIdx] / fluid_composition[activeCompIdx];
                        const ValueType rs = std::min(saturated_rs, max_possible_rs);
                        fluid_state.setRs(rs);
                    } else {
                        fluid_state.setRs(zero_value);
                    }
                }
                break;
            }
            case FluidSystem::gasPhaseIdx: {
                // TODO: we should use something else
                if constexpr (Indices::compositionSwitchIdx >= 0) {
                    if (both_oil_gas) {
                        const ValueType saturated_rv =
                            FluidSystem::saturatedVaporizationFactor(fluid_state, phaseIdx,
                                                                     fluid_state.pvtRegionIndex());
                        const unsigned oilCompIdx = FluidSystem::canonicalToActiveCompIdx(
                            FluidSystem::solventComponentIndex(FluidSystem::oilPhaseIdx));
                        const ValueType max_possible_rv =
                            fluid_composition[oilCompIdx] / fluid_composition[activeCompIdx];
                        const ValueType rv = std::min(saturated_rv, max_possible_rv);
                        fluid_state.setRv(rv);
                    } else {
                        fluid_state.setRv(zero_value);
                    }
                }
                if constexpr (Indices::waterSwitchIdx >= 0) {
                    if (both_water_gas && FluidSystem::enableVaporizedWater()) {
                        const ValueType saturated_rvw =
                            FluidSystem::saturatedVaporizationFactor(fluid_state, phaseIdx,
                                                                     fluid_state.pvtRegionIndex());
                        const unsigned waterCompIdx = FluidSystem::canonicalToActiveCompIdx(
                            FluidSystem::solventComponentIndex(FluidSystem::waterPhaseIdx));
                        const ValueType max_possible_rvw =
                            fluid_composition[waterCompIdx] / fluid_composition[activeCompIdx];
                        const ValueType rvw = std::min(saturated_rvw, max_possible_rvw);
                        fluid_state.setRvw(rvw);
                    } else {
                        fluid_state.setRvw(zero_value);
                    }
                }
                break;
            }
            case FluidSystem::waterPhaseIdx: {
                if constexpr (Indices::waterSwitchIdx >= 0) {
                    if (both_water_gas && FluidSystem::enableDissolvedGasInWater()) {
                        const ValueType saturated_rsw =
                            FluidSystem::saturatedDissolutionFactor(fluid_state, phaseIdx,
                                                                    fluid_state.pvtRegionIndex());
                        const unsigned gasCompIdx = FluidSystem::canonicalToActiveCompIdx(
                            FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                        const ValueType max_possible_rsw =
                            fluid_composition[gasCompIdx] / fluid_composition[activeCompIdx];
                        const ValueType rsw = std::min(saturated_rsw, max_possible_rsw);
                        fluid_state.setRsw(rsw);
                    } else {
                        fluid_state.setRsw(zero_value);
                    }
                }
                break;
            }
            default: {
                throw std::logic_error("Unhandled phase index " + std::to_string(phaseIdx));
            }
        }
        const auto& inv_b =
            FluidSystem::inverseFormationVolumeFactor(fluid_state, phaseIdx,
                                                      fluid_state.pvtRegionIndex());
        fluid_state.setInvB(phaseIdx, inv_b);
    }

    const bool both_water_gas_disgas =
        both_water_gas &&
        (FluidSystem::enableDissolvedGasInWater() || FluidSystem::enableVaporizedWater());

    std::vector<ValueType> saturation(FluidSystem::numPhases, zero_value);
    ValueType total_saturation {0.0};
    // calculate the saturation for all the phases
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }
        const bool water_gas_disgas_phase =
            both_water_gas_disgas &&
            (FluidSystem::waterPhaseIdx == phaseIdx || FluidSystem::gasPhaseIdx == phaseIdx);
        if (water_gas_disgas_phase) {
            // remove dissolved gas in water and vaporized water in gas
            const unsigned waterCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx);
            const unsigned gasCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx);
            // q_ws = q_wr * b_w + rvw * q_gr * b_g
            // q_gs = q_gr * b_g + rsw * q_wr * b_w
            // q_wr = 1 / (b_w * d) * (q_ws - rvw * q_gs)
            // q_gr = 1 / (b_g * d) * (q_gs - rsw * q_ws)
            // d = 1.0 - rsw * rvw
            const ValueType d = 1.0 - fluid_state.Rvw() * fluid_state.Rsw();
            if (d <= 0.0) {
                throw std::logic_error(
                    fmt::format("Problematic d value {} obtained for {}"
                                " during createBlackOilFluidState with rsw {}"
                                ", rvw {}.",
                                d, name, fluid_state.Rsw(), fluid_state.Rvw()));
            }
            if (FluidSystem::gasPhaseIdx == phaseIdx) {
                saturation[phaseIdx] =
                    (fluid_composition[gasCompIdx] -
                     fluid_state.Rsw() * fluid_composition[waterCompIdx]) /
                    (d * fluid_state.invB(phaseIdx));
            } else { // waterPhaseIdx
                saturation[phaseIdx] =
                    (fluid_composition[waterCompIdx] -
                     fluid_state.Rvw() * fluid_composition[gasCompIdx]) /
                    (d * fluid_state.invB(phaseIdx));
            }
            total_saturation += saturation[phaseIdx];
        } else if (!both_oil_gas || FluidSystem::waterPhaseIdx == phaseIdx) {
            const unsigned activeCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));
            saturation[phaseIdx] = fluid_composition[activeCompIdx] / fluid_state.invB(phaseIdx);
        } else {
            // remove dissolved gas and vaporized oil
            const unsigned oilCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx);
            const unsigned gasCompIdx =
                FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx);
            // q_os = q_or * b_o + rv * q_gr * b_g
            // q_gs = q_gr * b_g + rs * q_or * b_o
            // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)
            // d = 1.0 - rs * rv
            const ValueType d = 1.0 - fluid_state.Rv() * fluid_state.Rs();
            if (d <= 0.0) {
                throw std::logic_error(
                    fmt::format("Problematic d value {} obtained for {}"
                                " during createBlackOilFluidState with rs {}"
                                ", rv {}. Continue as if no dissolution (rs = 0) and"
                                " vaporization (rv = 0) for this connection.",
                                d, name, fluid_state.Rs(), fluid_state.Rv()));
            }
            if (FluidSystem::gasPhaseIdx == phaseIdx) {
                saturation[phaseIdx] =
                    (fluid_composition[gasCompIdx] -
                     fluid_state.Rs() * fluid_composition[oilCompIdx]) /
                    (d * fluid_state.invB(phaseIdx));
            } else if (FluidSystem::oilPhaseIdx == phaseIdx) {
                saturation[phaseIdx] =
                    (fluid_composition[oilCompIdx] -
                     fluid_state.Rv() * fluid_composition[gasCompIdx]) /
                    (d * fluid_state.invB(phaseIdx));
            }
            total_saturation += saturation[phaseIdx];
        }
    }

    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }
        fluid_state.setSaturation(phaseIdx, saturation[phaseIdx] / total_saturation);

        typename FluidSystem::template ParameterCache<ValueType> paramCache;
        paramCache.setRegionIndex(fluid_state.pvtRegionIndex());
        paramCache.updatePhase(fluid_state, phaseIdx);
        fluid_state.setDensity(phaseIdx, FluidSystem::density(fluid_state, paramCache, phaseIdx));
        if constexpr (has_energy) {
            fluid_state.setEnthalpy(phaseIdx,
                                    FluidSystem::enthalpy(fluid_state, paramCache, phaseIdx));
        }
    }

    return fluid_state;
}

} // namespace Opm

#endif // OPM_BLACKOIL_FLUID_STATE_UTILS_HPP
