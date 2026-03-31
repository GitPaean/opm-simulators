// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 SINTEF Digital

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
/*!
 * \file
 *
 * \copydoc Opm::FlowProblemComp
 */
#ifndef OPM_FLOW_PROBLEM_COMP_HPP
#define OPM_FLOW_PROBLEM_COMP_HPP


#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/FlowThresholdPressure.hpp>
#include <opm/simulators/flow/OutputCompositionalModule.hpp>
#include <opm/simulators/flow/equil/CompositionalEquil.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>

#include <opm/grid/utility/RegionMapping.hpp>

#include <algorithm>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup CompositionalSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class FlowProblemComp : public FlowProblem<TypeTag>
{
    // TODO: the naming of the Types will be adjusted
    using FlowProblemType = FlowProblem<TypeTag>;

    using typename FlowProblemType::Scalar;
    using typename FlowProblemType::Simulator;
    using typename FlowProblemType::GridView;
    using typename FlowProblemType::FluidSystem;
    using typename FlowProblemType::Vanguard;

    // might not be needed
    using FlowProblemType::dim;
    using FlowProblemType::dimWorld;

    using FlowProblemType::numPhases;
    using FlowProblemType::numComponents;

    using FlowProblemType::gasPhaseIdx;
    using FlowProblemType::oilPhaseIdx;
    using FlowProblemType::waterPhaseIdx;

    using typename FlowProblemType::Indices;
    using typename FlowProblemType::PrimaryVariables;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using typename FlowProblemType::Evaluation;
    using typename FlowProblemType::MaterialLaw;
    using typename FlowProblemType::RateVector;

    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using InitialFluidState = CompositionalFluidState<Scalar, FluidSystem>;
    using EclWriterType = EclWriter<TypeTag, OutputCompositionalModule<TypeTag> >;

public:
    using FlowProblemType::porosity;
    using FlowProblemType::pvtRegionIndex;

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        FlowProblemType::registerParameters();

        EclWriterType::registerParameters();

        // tighter tolerance is needed for compositional modeling here
        Parameters::SetDefault<Parameters::NewtonTolerance<Scalar>>(1e-7);
    }

    Opm::CompositionalConfig::EOSType getEosType() const
    {
        auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();
        return eclState.compositionalConfig().eosType(0);
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    explicit FlowProblemComp(Simulator& simulator)
        : FlowProblemType(simulator)
        , thresholdPressures_(simulator)
    {
        eclWriter_ = std::make_unique<EclWriterType>(simulator);
        enableEclOutput_ = Parameters::Get<Parameters::EnableEclOutput>();
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        // TODO: there should be room to remove duplication for this function,
        // but there is relatively complicated logic in the function calls in this function
        // some refactoring is needed for this function
        FlowProblemType::finishInit();

        auto& simulator = this->simulator();

        auto finishTransmissibilities = [updated = false, this]() mutable {
            if (updated) {
                return;
            }
            this->transmissibilities_.finishInit(
                [&vg = this->simulator().vanguard()](const unsigned int it) { return vg.gridIdxToEquilGridIdx(it); });
            updated = true;
        };
        // TODO: we might need to do the same with FlowProblemBlackoil for parallel

        finishTransmissibilities();

        if (enableEclOutput_) {
            eclWriter_->setTransmissibilities(&simulator.problem().eclTransmissibilities());
            std::function<unsigned int(unsigned int)> equilGridToGrid = [&simulator](unsigned int i) {
                return simulator.vanguard().gridEquilIdxToGridIdx(i);
            };
            eclWriter_->extractOutputTransAndNNC(equilGridToGrid);
        }

        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();

        // Set the start time of the simulation
        simulator.setStartTime(schedule.getStartTime());
        simulator.setEndTime(schedule.simTime(schedule.size() - 1));

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented by endEpisode(). The size of the initial time step and
        // length of the initial episode is set to zero for the same reason.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);

        // the "NOGRAV" keyword from Frontsim or setting the EnableGravity to false
        // disables gravity, else the standard value of the gravity constant at sea level
        // on earth is used
        this->gravity_ = 0.0;
        if (Parameters::Get<Parameters::EnableGravity>())
            this->gravity_[dim - 1] = 9.80665;
        if (!eclState.getInitConfig().hasGravity())
            this->gravity_[dim - 1] = 0.0;

        if (this->enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& tuning = schedule[0].tuning();
            this->initialTimeStepSize_ = tuning.TSINIT.has_value() ? tuning.TSINIT.value() : -1.0;
            this->maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        }

        this->initFluidSystem_();

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->maxOilSaturation_.resize(this->model().numGridDof(), 0.0);
        }

        this->readRockParameters_(simulator.vanguard().cellCenterDepths(), [&simulator](const unsigned idx) {
            std::array<int, dim> coords;
            simulator.vanguard().cartesianCoordinate(idx, coords);
            std::ranges::transform(coords, coords.begin(),
                                   [](const auto c) { return c + 1; });
            return coords;
        });
        FlowProblemType::readMaterialParameters_();
        FlowProblemType::readThermalParameters_();

        // write the static output files (EGRID, INIT)
        if (enableEclOutput_) {
            eclWriter_->writeInit();
        }

        const auto& initconfig = eclState.getInitConfig();
        if (initconfig.restartRequested())
            readEclRestartSolution_();
        else
            this->readInitialCondition_();

        FlowProblemType::updatePffDofData_();

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            this->polymer_.maxAdsorption.resize(numElements, 0.0);
        }

        /* readBoundaryConditions_();

        // compute and set eq weights based on initial b values
        computeAndSetEqWeights_();

        if (enableDriftCompensation_) {
            drift_.resize(this->model().numGridDof());
            drift_ = 0.0;
        } */

        // TODO: check wether the following can work with compostional
        if (this->enableVtkOutput_() && eclState.getIOConfig().initOnly()) {
            simulator.setTimeStepSize(0.0);
            FlowProblemType::writeOutput(true);
        }

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested()) {
            simulator.startNextEpisode(schedule.seconds(1));
            simulator.setEpisodeIndex(0);
            simulator.setTimeStepIndex(0);
        }
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep() override
    {
        FlowProblemType::endTimeStep();

        // after the solution is updated, the values in output module also needs to be updated
        this->eclWriter_->mutableOutputModule().invalidateLocalData();

        // For CpGrid with LGRs, ecl/vtk output is not supported yet.
        const auto& grid = this->simulator().vanguard().gridView().grid();

        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;
        if (!isCpGrid || (grid.maxLevel() == 0)) {
            this->eclWriter_->evalSummaryState(! this->episodeWillBeOver());
        }
    }

    void writeReports(const SimulatorTimer& timer) {
        if (enableEclOutput_){
            eclWriter_->writeReports(timer);
        }
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(bool verbose) override
    {
        FlowProblemType::writeOutput(verbose);

        if (! this->enableEclOutput_) {
            return;
        }

        const auto isSubStep = !this->episodeWillBeOver();

        if (!isSubStep || Parameters::Get<Parameters::EnableWriteAllSolutions>()) {
            auto localCellData = data::Solution {};

            this->eclWriter_->writeOutput(std::move(localCellData), isSubStep);
        }
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * Reservoir simulation uses no-flow conditions as default for all boundaries.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned /* timeIdx */) const
    {
        OPM_TIMEBLOCK_LOCAL(eclProblemBoundary, Subsystem::Assembly);
        if (!context.intersection(spaceIdx).boundary())
            return;

        values.setNoFlow();

        if (this->nonTrivialBoundaryConditions()) {
            throw std::logic_error("boundary condition is not supported by compostional modeling yet");
        }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        const auto& initial_fs = initialFluidStates_[globalDofIdx];
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        for (unsigned p = 0; p < numPhases; ++p) { // TODO: assuming the phaseidx continuous
            // pressure
            fs.setPressure(p, initial_fs.pressure(p));

            // saturation
            fs.setSaturation(p, initial_fs.saturation(p));

            // temperature
            fs.setTemperature(initial_fs.temperature(p));
        }


        if (!zmf_initialization_) {
            for (unsigned p = 0; p < numPhases; ++p) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    fs.setMoleFraction(p, compIdx, initial_fs.moleFraction(p, compIdx));
                }
            }

            {
                const auto& eos_type = getEosType();
                typename FluidSystem::template ParameterCache<Scalar> paramCache(eos_type);
                paramCache.updatePhase(fs, FluidSystem::oilPhaseIdx);
                paramCache.updatePhase(fs, FluidSystem::gasPhaseIdx);
                fs.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::oilPhaseIdx));
                fs.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::gasPhaseIdx));
            }
            // determine the component fractions
            Dune::FieldVector<Scalar, numComponents> z(0.0);
            Scalar sumMoles = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (Indices::waterEnabled && phaseIdx == static_cast<unsigned int>(waterPhaseIdx)){
                    continue;
                }
                const auto saturation = fs.saturation(phaseIdx);
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar tmp = fs.molarity(phaseIdx, compIdx) * saturation;
                    tmp = max(tmp, 1e-8);
                    z[compIdx] += tmp;
                    sumMoles += tmp;
                }
            }
            z /= sumMoles;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setMoleFraction(compIdx, z[compIdx]);
            }
        } else {
            // TODO: should we normalize the input?
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                fs.setMoleFraction(compIdx, initial_fs.moleFraction(compIdx));
            }
        }

        // Set initial K and L
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const auto& Ktmp = fs.wilsonK_(compIdx);
            fs.setKvalue(compIdx, Ktmp);
        }

        const Scalar& Ltmp = -1.0;
        fs.setLvalue(Ltmp);

        values.assignNaive(fs);
    }

    void addToSourceDense(RateVector&, unsigned, unsigned) const override
    {
        // we do nothing for now
    }

    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    std::vector<InitialFluidState>& initialFluidStates()
    { return initialFluidStates_; }

    const std::vector<InitialFluidState>& initialFluidStates() const
    { return initialFluidStates_; }

    const FlowThresholdPressure<TypeTag>& thresholdPressure() const
    {
        assert( !thresholdPressures_.enableThresholdPressure() &&
                " Threshold Pressures are not supported by compostional simulation ");
        return thresholdPressures_;
    }

    const EclWriterType& eclWriter() const
    { return *eclWriter_; }

    EclWriterType& eclWriter()
    { return *eclWriter_; }

    // TODO: do we need this one?
    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<FlowProblemType&>(*this));
        serializer(*eclWriter_);
    }
protected:

    void updateExplicitQuantities_(int /* episodeIdx*/, int /* timeStepSize */, bool /* first_step_after_restart */) override
    {
        // we do nothing here for now
    }

    void readEquilInitialCondition_() override
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& initConfig = eclState.getInitConfig();

        if (!initConfig.hasEquil()) {
            throw std::domain_error(
                "Compositional equilibration requires EQUIL keyword in the deck");
        }

        const auto& tables = eclState.getTableManager();
        const auto& equil = initConfig.getEquil();
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
        const std::size_t numDof = this->model().numGridDof();

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
                    // No ZMFVD for this region — use uniform composition from first region's shallowest entry
                    throw std::runtime_error(
                        "ZMFVD table not provided for EQUIL region "
                        + std::to_string(r + 1)
                        + ". Compositional equilibration requires ZMFVD for all regions.");
                }
            }
        } else {
            // No ZMFVD keyword at all — check if ZMF field property is available
            if (eclState.fieldProps().has_double("ZMF")) {
                // Use a uniform composition from the first cell as fallback
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

        const int numPressurePoints = this->numPressurePointsEquil();
        const Scalar grav = this->gravity_[dim - 1];

        // Prepare output storage
        initialFluidStates_.resize(numDof);
        zmf_initialization_ = true;  // We store overall compositions (z_i)

        // We also need per-cell pressure storage like the Black-Oil code
        // For compositional, pressures and saturations go into initialFluidStates_ directly.

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

            // Build the compositional pressure table
            // We need an EquilReg-like object for datum/contacts.
            // For compositional, we only need datum depth, datum pressure, GOC, and Pcgo.
            // We create a minimal EquilReg.
            const auto rsFunc = std::make_shared<EQUIL::Miscibility::NoMixing<Scalar>>();
            const auto rvFunc = std::make_shared<EQUIL::Miscibility::NoMixing<Scalar>>();
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
                const Scalar temp = cellTemperature[cell];

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
                    // Water pressure approximation: P_wat ≈ P_oil at WOC, then hydrostatic.
                    // For now use oil pressure as proxy (TODO: proper water pressure table)
                    fs.setPressure(waterPhaseIdx, pOil);
                }

                // Set saturations
                fs.setSaturation(oilPhaseIdx, result.So);
                fs.setSaturation(gasPhaseIdx, result.Sg);
                if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    // TODO: proper water saturation from water pressure table + Pc inversion
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

    void readEclRestartSolution_()
    {
        throw std::logic_error("Restarting is not supported by compositional modeling yet");
    }

    void readExplicitInitialCondition_() override
    {
        readExplicitInitialConditionCompositional_();
    }

    void readExplicitInitialConditionCompositional_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& fp = eclState.fieldProps();
        const bool has_pressure = fp.has_double("PRESSURE");
        if (!has_pressure)
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                     "keyword if the model is initialized explicitly");

        const bool has_xmf = fp.has_double("XMF");
        const bool has_ymf = fp.has_double("YMF");
        const bool has_zmf = fp.has_double("ZMF");
        if ( !has_zmf && !(has_xmf && has_ymf) ) {
            throw std::runtime_error("The ECL input file requires the presence of ZMF or XMF and YMF "
                                     "keyword if the model is initialized explicitly");
        }

        if (has_zmf && (has_xmf || has_ymf)) {
            throw std::runtime_error("The ECL input file can not handle explicit initialization "
                                     "with both ZMF and XMF or YMF");
        }

        if (has_xmf != has_ymf) {
            throw std::runtime_error("The ECL input file needs XMF and YMF combined to do the explicit "
                                     "initializtion when using XMF or YMF");
        }

        const bool has_temp = fp.has_double("TEMPI");

        // const bool has_gas = fp.has_double("SGAS");
        assert(fp.has_double("SGAS"));

        std::size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> soilData;
        std::vector<double> pressureData;
        std::vector<double> tempiData;

        const bool water_active = FluidSystem::phaseIsActive(waterPhaseIdx);
        const bool gas_active = FluidSystem::phaseIsActive(gasPhaseIdx);
        const bool oil_active = FluidSystem::phaseIsActive(oilPhaseIdx);

        if (water_active && Indices::numPhases > 2)
            waterSaturationData = fp.get_double("SWAT");
        else
            waterSaturationData.resize(numDof);

        pressureData = fp.get_double("PRESSURE");

        if (has_temp) {
            tempiData = fp.get_double("TEMPI");
        } else {
            ; // TODO: throw?
        }

        if (gas_active) // && FluidSystem::phaseIsActive(oilPhaseIdx))
            gasSaturationData = fp.get_double("SGAS");
        else
            gasSaturationData.resize(numDof);

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];
            // dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            Scalar temperatureLoc = tempiData[dofIdx];
            assert(std::isfinite(temperatureLoc) && temperatureLoc > 0);
            dofFluidState.setTemperature(temperatureLoc);

            if (gas_active) {
                dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                            gasSaturationData[dofIdx]);
            }
            if (oil_active) {
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);
            }
            if (water_active) {
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);
            }

            //////
            // set phase pressures
            //////
            const Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // TODO: zero capillary pressure for now
            const std::array<Scalar, numPhases> pc = {0};
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                if (Indices::oilEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                else if (Indices::gasEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                else if (Indices::waterEnabled)
                    // single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
            }

            if (has_xmf && has_ymf) {
                const auto& xmfData = fp.get_double("XMF");
                const auto& ymfData = fp.get_double("YMF");
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const std::size_t data_idx = compIdx * numDof + dofIdx;
                    const Scalar xmf = xmfData[data_idx];
                    const Scalar ymf = ymfData[data_idx];

                    dofFluidState.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, xmf);
                    dofFluidState.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, ymf);
                }
            }

            if (has_zmf) {
                zmf_initialization_ = true;
                const auto& zmfData = fp.get_double("ZMF");
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const std::size_t data_idx = compIdx * numDof + dofIdx;
                    const Scalar zmf = zmfData[data_idx];
                    dofFluidState.setMoleFraction(compIdx, zmf);

                    if (gas_active) {
                        const auto ymf = (dofFluidState.saturation(FluidSystem::gasPhaseIdx) > 0.) ? zmf : Scalar{0};
                        dofFluidState.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, ymf);
                    }
                    if (oil_active) {
                        const auto xmf = (dofFluidState.saturation(FluidSystem::oilPhaseIdx) > 0.) ? zmf : Scalar{0};
                        dofFluidState.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, xmf);
                    }
                }
            }
        }
    }

private:

    void handleSolventBC(const BCProp::BCFace& /* bc */, RateVector& /* rate */) const override
    {
        throw std::logic_error("solvent is disabled for compositional modeling and you're trying to add solvent to BC");
    }

    void handlePolymerBC(const BCProp::BCFace& /* bc */, RateVector& /* rate */) const override
    {
        throw std::logic_error("polymer is disabled for compositional modeling and you're trying to add polymer to BC");
    }

    void handleMicrBC(const BCProp::BCFace& /* bc */, RateVector& /* rate */) const override
    {
        throw std::logic_error("MICP is disabled for compositional modeling and you're trying to add microbes to BC");
    }

    void handleOxygBC(const BCProp::BCFace& /* bc */, RateVector& /* rate */) const override
    {
        throw std::logic_error("MICP is disabled for compositional modeling and you're trying to add oxygen to BC");
    }

    void handleUreaBC(const BCProp::BCFace& /* bc */, RateVector& /* rate */) const override
    {
        throw std::logic_error("MICP is disabled for compositional modeling and you're trying to add urea to BC");
    }

    FlowThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<InitialFluidState> initialFluidStates_;

    bool zmf_initialization_ {false};

    bool enableEclOutput_{false};
    std::unique_ptr<EclWriterType> eclWriter_;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_COMP_HPP
