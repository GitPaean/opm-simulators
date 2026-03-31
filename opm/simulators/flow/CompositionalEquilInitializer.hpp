// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \copydoc Opm::CompositionalEquilInitializer
 */
#ifndef OPM_COMPOSITIONAL_EQUIL_INITIALIZER_HPP
#define OPM_COMPOSITIONAL_EQUIL_INITIALIZER_HPP

#include <opm/simulators/flow/equil/CompositionalEquil.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup CompositionalSimulator
 *
 * \brief Computes the initial condition for compositional models based
 *        on the EQUIL keyword from ECL.
 *
 * Uses flash-based equilibration: pressure-vs-depth tables are built
 * by integrating dP/dz = rho(z,P)*g where the density comes from an
 * EoS flash at each (depth, pressure, composition) sample.  Per-cell
 * saturations and phase compositions are then determined by a final
 * flash at each cell centre.
 *
 * Analogous to EquilInitializer for Black-Oil models.
 */
template <class TypeTag>
class CompositionalEquilInitializer
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr int numPhases     = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int oilPhaseIdx   = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx   = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr int dimWorld      = GridView::dimensionworld;

    using InitialFluidState = CompositionalFluidState<Scalar, FluidSystem>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using EOSType = CompositionalConfig::EOSType;

public:
    /*!
     * \brief Construct the initializer and compute the initial state.
     *
     * \param simulator  The simulator object (provides access to grid,
     *                   eclState, gravity, etc.).
     */
    explicit CompositionalEquilInitializer(const Simulator& simulator)
        : simulator_(simulator)
    {
        const auto& vanguard  = simulator.vanguard();
        const auto& eclState  = vanguard.eclState();
        const auto& initConfig = eclState.getInitConfig();

        if (!initConfig.hasEquil()) {
            throw std::domain_error(
                "Compositional equilibration requires EQUIL keyword in the deck");
        }

        const auto& tables = eclState.getTableManager();
        const auto& equil  = initConfig.getEquil();
        const std::vector<EquilRecord> rec(equil.begin(), equil.end());
        const std::size_t numRegions = rec.size();

        // --- EQLNUM region mapping ---
        const auto& gridView = vanguard.gridView();
        std::vector<int> eqlnum(gridView.size(/*codim=*/0), 0);
        if (eclState.fieldProps().has_int("EQLNUM")) {
            const auto& e = eclState.fieldProps().get_int("EQLNUM");
            std::transform(e.begin(), e.end(), eqlnum.begin(),
                           [](int n) { return n - 1; });
        }

        // --- EoS type ---
        const auto eosType = eclState.compositionalConfig().eosType(0);

        // --- Cell center depths ---
        const auto& cellCenterDepth = vanguard.cellCenterDepths();
        const std::size_t numDof = gridView.size(/*codim=*/0);

        // --- Temperature vs depth tables ---
        std::vector<TabulatedFunction> tempVdTable(numRegions);
        if (!tables.hasTables("RTEMPVD")) {
            const Scalar rtemp = static_cast<Scalar>(tables.rtemp());
            std::vector<Scalar> x = {0.0, 1.0};
            std::vector<Scalar> y = {rtemp, rtemp};
            for (auto& t : tempVdTable) {
                t.setXYContainers(x, y);
            }
        } else {
            const auto& tempvdTables = tables.getRtempvdTables();
            for (std::size_t i = 0; i < tempvdTables.size(); ++i) {
                const auto& tvd = tempvdTables.template getTable<RtempvdTable>(i);
                tempVdTable[i].setXYContainers(tvd.getDepthColumn(),
                                               tvd.getTemperatureColumn());
            }
        }

        // --- Composition vs depth (ZMFVD) ---
        using CompositionVsDepth = EQUIL::Details::CompositionVsDepth<Scalar>;
        const auto& zmfvdTables = tables.getZmfvdTables();
        const int numComp = static_cast<int>(numComponents);

        std::vector<CompositionVsDepth> compositionVsDepth;
        compositionVsDepth.reserve(numRegions);

        if (!zmfvdTables.empty()) {
            for (std::size_t r = 0; r < numRegions; ++r) {
                if (r < zmfvdTables.size()) {
                    const auto& zmfvd = zmfvdTables.template getTable<ZmfvdTable>(r);
                    const auto& depthCol = zmfvd.getDepthColumn();

                    std::vector<Scalar> depths(depthCol.size());
                    for (std::size_t j = 0; j < depthCol.size(); ++j) {
                        depths[j] = static_cast<Scalar>(depthCol[j]);
                    }

                    std::vector<std::vector<Scalar>> compositions(numComp);
                    for (int c = 0; c < numComp; ++c) {
                        const auto& col = zmfvd.getMoleFractionColumn(c);
                        compositions[c].resize(col.size());
                        for (std::size_t j = 0; j < col.size(); ++j) {
                            compositions[c][j] = static_cast<Scalar>(col[j]);
                        }
                    }
                    compositionVsDepth.emplace_back(depths, compositions, numComp);
                } else {
                    throw std::runtime_error(
                        "ZMFVD table not provided for EQUIL region "
                        + std::to_string(r + 1)
                        + ". Compositional equilibration requires ZMFVD for all regions.");
                }
            }
        } else {
            // No ZMFVD keyword at all — check if ZMF field property is available
            if (eclState.fieldProps().has_double("ZMF")) {
                const auto& zmfData = eclState.fieldProps().get_double("ZMF");
                std::vector<Scalar> zi(numComp, Scalar{0});
                for (int c = 0; c < numComp; ++c) {
                    zi[c] = static_cast<Scalar>(zmfData[c * numDof]);
                }
                for (std::size_t r = 0; r < numRegions; ++r) {
                    compositionVsDepth.emplace_back(zi);
                }
                OpmLog::warning("No ZMFVD tables found. Using ZMF from first cell as "
                                "constant composition for equilibration.");
            } else {
                throw std::runtime_error(
                    "Compositional equilibration requires ZMFVD keyword or ZMF field data");
            }
        }

        // --- Compute cell temperatures (for later use) ---
        std::vector<Scalar> cellTemperature(numDof);
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            const int reg = eqlnum[dofIdx];
            cellTemperature[dofIdx] = tempVdTable[reg].eval(
                cellCenterDepth[dofIdx], /*extrapolate=*/true);
        }

        // --- Build pressure tables and compute initial state per region ---
        using CompPressureTable =
            EQUIL::Details::CompositionalPressureTable<FluidSystem, EQUIL::EquilReg<Scalar>>;
        using CompPhaseSat =
            EQUIL::Details::CompositionalPhaseSaturations<FluidSystem>;

        const int numPressurePoints = simulator.problem().numPressurePointsEquil();
        const Scalar grav = simulator.problem().gravity()[dimWorld - 1];

        // Prepare output storage
        initialFluidStates_.resize(numDof);

        for (std::size_t r = 0; r < numRegions; ++r) {
            // Collect cells in this region
            std::vector<std::size_t> regionCells;
            for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                if (static_cast<std::size_t>(eqlnum[dofIdx]) == r) {
                    regionCells.push_back(dofIdx);
                }
            }

            if (regionCells.empty()) {
                continue;
            }

            // Compute vertical extent of this region
            Scalar zMin = std::numeric_limits<Scalar>::max();
            Scalar zMax = std::numeric_limits<Scalar>::lowest();
            for (const auto& cell : regionCells) {
                const Scalar depth = cellCenterDepth[cell];
                zMin = std::min(zMin, depth);
                zMax = std::max(zMax, depth);
            }

            // Extend span to include contacts — but only for active phases.
            // For a two-phase oil+gas compositional system, the WOC is
            // irrelevant and should NOT extend the span.  Each depth in
            // the span triggers a full flash evaluation in the RK4 solver,
            // so keeping the span tight is important for performance.
            const Scalar zgoc = rec[r].gasOilContactDepth();
            const Scalar zwoc = rec[r].waterOilContactDepth();

            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                zMin = std::min(zMin, zgoc);
                zMax = std::max(zMax, zgoc);
            }
            if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                zMin = std::min(zMin, zwoc);
                zMax = std::max(zMax, zwoc);
            }

            // Add a small buffer
            const Scalar dz = (zMax - zMin) * 0.01;
            std::array<Scalar, 2> vspan = {{ zMin - dz, zMax + dz }};

            // Build EquilReg for this region (compositional uses NoMixing)
            const auto rsFunc  = std::make_shared<EQUIL::Miscibility::NoMixing<Scalar>>();
            const auto rvFunc  = std::make_shared<EQUIL::Miscibility::NoMixing<Scalar>>();
            const auto rvwFunc = std::make_shared<EQUIL::Miscibility::NoMixing<Scalar>>();
            TabulatedFunction saltVdTable;
            {
                std::vector<Scalar> x = {0.0, 1.0};
                std::vector<Scalar> y = {0.0, 0.0};
                saltVdTable.setXYContainers(x, y);
            }

            const auto eqreg = EQUIL::EquilReg<Scalar>{
                rec[r], rsFunc, rvFunc, rvwFunc,
                tempVdTable[r], saltVdTable, /*pvtIdx=*/0
            };

            // Create and equilibrate the compositional pressure table
            auto compPtable = CompPressureTable{
                grav, numPressurePoints, compositionVsDepth[r], eosType
            };
            compPtable.equilibrate(eqreg, vspan);

            // Create the phase saturation calculator
            auto compPhaseSat = CompPhaseSat{eosType};

            // For each cell in this region, compute initial state
            for (const auto& cell : regionCells) {
                auto& fs = initialFluidStates_[cell];

                const Scalar depth = cellCenterDepth[cell];
                const Scalar temp  = cellTemperature[cell];

                // Look up pressures from the pre-computed tables
                const Scalar pOil = compPtable.oil(depth);
                const Scalar pGas = compPtable.gas(depth);

                // Look up composition at this depth
                const auto zi = compositionVsDepth[r].eval(depth);

                // Flash to get saturations and phase compositions
                const auto result = compPhaseSat.compute(pOil, pGas, temp, zi);

                // Set temperature
                fs.setTemperature(temp);

                // Set pressures
                fs.setPressure(oilPhaseIdx, pOil);
                fs.setPressure(gasPhaseIdx, pGas);
                if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    // Water pressure approximation: P_wat ≈ P_oil
                    // TODO: proper water pressure table
                    fs.setPressure(waterPhaseIdx, pOil);
                }

                // Set saturations
                fs.setSaturation(oilPhaseIdx, result.So);
                fs.setSaturation(gasPhaseIdx, result.Sg);
                if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    // TODO: proper water saturation from Pc inversion
                    fs.setSaturation(waterPhaseIdx, Scalar{0});
                }

                // Set overall mole fractions (z_i)
                for (unsigned c = 0; c < numComponents; ++c) {
                    fs.setMoleFraction(c, zi[c]);
                }

                // Set phase compositions from flash result
                for (unsigned c = 0; c < numComponents; ++c) {
                    fs.setMoleFraction(oilPhaseIdx, c, result.x[c]);
                    fs.setMoleFraction(gasPhaseIdx, c, result.y[c]);
                }

                // Set K-values and L-value from flash
                fs.setLvalue(result.L);
                for (unsigned c = 0; c < numComponents; ++c) {
                    if (result.x[c] > 1e-20) {
                        fs.setKvalue(c, result.y[c] / result.x[c]);
                    } else {
                        fs.setKvalue(c, fs.wilsonK_(c));
                    }
                }

                // Compute and set densities using the EoS
                {
                    typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType);
                    if (result.So > 0) {
                        paramCache.updatePhase(fs, oilPhaseIdx);
                        fs.setDensity(oilPhaseIdx,
                                      FluidSystem::density(fs, paramCache, oilPhaseIdx));
                    }
                    if (result.Sg > 0) {
                        paramCache.updatePhase(fs, gasPhaseIdx);
                        fs.setDensity(gasPhaseIdx,
                                      FluidSystem::density(fs, paramCache, gasPhaseIdx));
                    }
                }
            }
        }

        OpmLog::info("Compositional equilibration completed successfully");
    }

    /*!
     * \brief Return the initial thermodynamic state which should be used
     *        as the initial condition.
     *
     * This is supposed to correspond to hydrostatic conditions.
     */
    const InitialFluidState& initialFluidState(unsigned elemIdx) const
    {
        return initialFluidStates_[elemIdx];
    }

private:
    const Simulator& simulator_;
    std::vector<InitialFluidState> initialFluidStates_;
};

} // namespace Opm

#endif // OPM_COMPOSITIONAL_EQUIL_INITIALIZER_HPP
