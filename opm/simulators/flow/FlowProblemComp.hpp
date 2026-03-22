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

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

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
    using typename FlowProblemType::TabulatedFunction;

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
        const auto& tables = eclState.getTableManager();
        const auto eos_type = getEosType();
        const Scalar grav = this->gravity_[dim - 1];
        const auto& cellCenterDepths = vanguard.cellCenterDepths();
        const std::size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        const bool water_active = FluidSystem::phaseIsActive(waterPhaseIdx);
        const bool gas_active = FluidSystem::phaseIsActive(gasPhaseIdx);
        const bool oil_active = FluidSystem::phaseIsActive(oilPhaseIdx);

        // --- Temperature vs depth ---
        const Scalar refTemp = tables.rtemp();
        TabulatedFunction tempVdFunc;
        if (tables.hasTables("RTEMPVD")) {
            const auto& tempvdTables = tables.getRtempvdTables();
            const auto& rtempvdTable = tempvdTables.template getTable<RtempvdTable>(0);
            tempVdFunc.setXYContainers(rtempvdTable.getDepthColumn(),
                                       rtempvdTable.getTemperatureColumn());
        } else {
            std::vector<Scalar> defaultDepths = {Scalar(0), Scalar(1e6)};
            std::vector<Scalar> defaultTemps = {refTemp, refTemp};
            tempVdFunc.setXYContainers(defaultDepths, defaultTemps);
        }

        // --- Overall composition vs depth ---
        // If ZMFVD is present, create per-component interpolation tables;
        // otherwise use uniform molar fractions as a default.
        const bool hasZmfvd = tables.hasTables("ZMFVD");
        std::vector<TabulatedFunction> compVdFuncs(numComponents);

        if (hasZmfvd) {
            const auto& zmfvdTables = tables.getZmfvdTables();
            const auto& zmfvdTable = zmfvdTables.template getTable<ZmfvdTable>(0);
            const auto& depthCol = zmfvdTable.getDepthColumn();
            for (unsigned c = 0; c < numComponents; ++c) {
                const auto& fracCol = zmfvdTable.getMoleFractionColumn(c);
                compVdFuncs[c].setXYContainers(depthCol, fracCol);
            }
            OpmLog::info("Compositional equilibration: using ZMFVD for "
                         "depth-dependent composition.");
        } else {
            const Scalar uniformFrac = Scalar(1.0) / static_cast<Scalar>(numComponents);
            std::vector<Scalar> defaultDepths = {Scalar(0), Scalar(1e6)};
            for (unsigned c = 0; c < numComponents; ++c) {
                std::vector<Scalar> defaultFracs = {uniformFrac, uniformFrac};
                compVdFuncs[c].setXYContainers(defaultDepths, defaultFracs);
            }
            OpmLog::info("Compositional equilibration: using uniform molar "
                         "composition for pressure profile computation.");
        }

        // Helper: return overall composition at a given depth
        auto compositionAtDepth = [&](Scalar depth, std::vector<Scalar>& comp) {
            Scalar sum = 0;
            for (unsigned c = 0; c < numComponents; ++c) {
                comp[c] = std::max(Scalar(0), compVdFuncs[c].eval(depth, true));
                sum += comp[c];
            }
            // Normalize so that mole fractions sum to 1
            if (sum > Scalar(0)) {
                for (unsigned c = 0; c < numComponents; ++c)
                    comp[c] /= sum;
            }
        };

        // --- EQUIL data ---
        const auto& initConfig = eclState.getInitConfig();
        const auto& equil = initConfig.getEquil();
        // Process first equilibration region
        // TODO: support multiple equilibration regions
        const auto& rec = equil.getRecord(0);
        const Scalar datumDepth = rec.datumDepth();
        const Scalar datumPressure = rec.datumDepthPressure();

        // Contact depths: use extreme values when phase pair is not active
        constexpr Scalar noContactBelow = Scalar(1e20);
        constexpr Scalar noContactAbove = Scalar(-1e20);
        const Scalar wocDepth = (water_active && oil_active)
            ? rec.waterOilContactDepth() : noContactBelow;
        const Scalar gocDepth = (gas_active && oil_active)
            ? rec.gasOilContactDepth() : noContactAbove;
        const Scalar wocPc = (water_active && oil_active)
            ? rec.waterOilContactCapillaryPressure() : Scalar(0);
        const Scalar gocPc = (gas_active && oil_active)
            ? rec.gasOilContactCapillaryPressure() : Scalar(0);

        // --- Depth range ---
        Scalar minDepth = cellCenterDepths[0];
        Scalar maxDepth = cellCenterDepths[0];
        for (std::size_t i = 1; i < numDof; ++i) {
            minDepth = std::min(minDepth, cellCenterDepths[i]);
            maxDepth = std::max(maxDepth, cellCenterDepths[i]);
        }
        if (gas_active && oil_active)
            minDepth = std::min(minDepth, gocDepth);
        if (water_active && oil_active)
            maxDepth = std::max(maxDepth, wocDepth);
        minDepth = std::min(minDepth, datumDepth) - Scalar(1.0);
        maxDepth = std::max(maxDepth, datumDepth) + Scalar(1.0);
        const Scalar depthRange = maxDepth - minDepth;

        // --- Density computation ---
        // For oil/gas: uses cubic EOS with overall composition.
        // The ParameterCache selects the liquid root for oil and the gas root
        // for gas, so no flash is needed for the pressure integration.
        auto hydrocarbonDensity = [&](Scalar depth, Scalar pressure,
                                      unsigned phaseIdx) -> Scalar {
            CompositionalFluidState<Scalar, FluidSystem> fs;
            fs.setTemperature(tempVdFunc.eval(depth, true));
            for (unsigned p = 0; p < numPhases; ++p)
                fs.setPressure(p, pressure);
            // Evaluate composition at the given depth
            std::vector<Scalar> comp(numComponents);
            compositionAtDepth(depth, comp);
            // Set composition for all hydrocarbon phases (ParameterCache needs both)
            for (unsigned c = 0; c < numComponents; ++c) {
                if (oil_active)
                    fs.setMoleFraction(oilPhaseIdx, c, comp[c]);
                if (gas_active)
                    fs.setMoleFraction(gasPhaseIdx, c, comp[c]);
            }
            typename FluidSystem::template ParameterCache<Scalar> paramCache(eos_type);
            paramCache.updatePhase(fs, phaseIdx);
            return FluidSystem::density(fs, paramCache, phaseIdx);
        };

        // For water: uses water PVT tables (no EOS needed)
        auto waterDensityFn = [&](Scalar depth, Scalar pressure) -> Scalar {
            CompositionalFluidState<Scalar, FluidSystem> fs;
            fs.setTemperature(tempVdFunc.eval(depth, true));
            fs.setPressure(waterPhaseIdx, pressure);
            typename FluidSystem::template ParameterCache<Scalar> paramCache(eos_type);
            return FluidSystem::density(fs, paramCache, waterPhaseIdx);
        };

        // --- RK4 pressure integration: dP/dz = rho(z,P) * g ---
        const int numPressureSamples = 2000;
        constexpr Scalar minRangeLength = Scalar(1e-6);

        auto integratePhase = [&](Scalar startDepth, Scalar startPressure,
                                  Scalar endDepth, auto densityFn)
        {
            const Scalar rangeLen = std::abs(endDepth - startDepth);
            if (rangeLen < minRangeLength) {
                return std::make_pair(
                    std::vector<Scalar>{startDepth},
                    std::vector<Scalar>{startPressure});
            }
            const int nSteps = std::max(10,
                static_cast<int>(rangeLen / (depthRange + Scalar(1)) * numPressureSamples));
            const Scalar dz = (endDepth - startDepth) / nSteps;

            std::vector<Scalar> depths(nSteps + 1);
            std::vector<Scalar> pressures(nSteps + 1);
            depths[0] = startDepth;
            pressures[0] = startPressure;

            for (int i = 0; i < nSteps; ++i) {
                const Scalar z = depths[i];
                const Scalar P = pressures[i];
                const Scalar k1 = densityFn(z, P) * grav * dz;
                const Scalar k2 = densityFn(z + Scalar(0.5) * dz,
                                            P + Scalar(0.5) * k1) * grav * dz;
                const Scalar k3 = densityFn(z + Scalar(0.5) * dz,
                                            P + Scalar(0.5) * k2) * grav * dz;
                const Scalar k4 = densityFn(z + dz, P + k3) * grav * dz;
                depths[i + 1] = z + dz;
                pressures[i + 1] = P + (k1 + Scalar(2) * k2
                                        + Scalar(2) * k3 + k4) / Scalar(6);
            }
            return std::make_pair(std::move(depths), std::move(pressures));
        };

        // Builds a pressure lookup table by integrating up and down from a
        // known (depth, pressure) pair and merging the results.
        auto buildPressureTable = [&](Scalar startDepth, Scalar startPressure,
                                      auto densityFn) -> TabulatedFunction
        {
            auto [dUp, pUp] = integratePhase(startDepth, startPressure,
                                             minDepth, densityFn);
            auto [dDown, pDown] = integratePhase(startDepth, startPressure,
                                                 maxDepth, densityFn);

            // Merge: reversed upward (skip first = startDepth) + downward
            std::vector<Scalar> allDepths;
            std::vector<Scalar> allPressures;
            allDepths.reserve(dUp.size() + dDown.size());
            allPressures.reserve(pUp.size() + pDown.size());

            for (int i = static_cast<int>(dUp.size()) - 1; i >= 1; --i) {
                allDepths.push_back(dUp[i]);
                allPressures.push_back(pUp[i]);
            }
            for (std::size_t i = 0; i < dDown.size(); ++i) {
                allDepths.push_back(dDown[i]);
                allPressures.push_back(pDown[i]);
            }

            // Ensure at least 2 points for interpolation
            if (allDepths.size() < 2) {
                allDepths = {minDepth, maxDepth};
                allPressures = {startPressure, startPressure};
            }

            TabulatedFunction func;
            func.setXYContainers(allDepths, allPressures);
            return func;
        };

        // --- Build phase pressure tables ---
        // Strategy depends on datum depth relative to contacts.
        // The convention is: Pcow = Po - Pw, Pcgo = Pg - Po
        TabulatedFunction oilPressFunc, gasPressFunc, waterPressFunc;
        bool oilBuilt = false, gasBuilt = false, waterBuilt = false;

        const bool datumInGas = oil_active && gas_active && (datumDepth < gocDepth);
        const bool datumInWater = oil_active && water_active && (datumDepth > wocDepth);

        if (oil_active) {
            auto oilDensityFn = [&](Scalar depth, Scalar press) {
                return hydrocarbonDensity(depth, press, oilPhaseIdx);
            };

            if (datumInGas) {
                // Datum in gas zone: build gas first, derive oil from GOC
                auto gasDensityFn = [&](Scalar depth, Scalar press) {
                    return hydrocarbonDensity(depth, press, gasPhaseIdx);
                };
                gasPressFunc = buildPressureTable(datumDepth, datumPressure,
                                                  gasDensityFn);
                gasBuilt = true;
                const Scalar oilPAtGoc = gasPressFunc.eval(gocDepth, true) - gocPc;
                oilPressFunc = buildPressureTable(gocDepth, oilPAtGoc, oilDensityFn);
                oilBuilt = true;
            } else if (datumInWater) {
                // Datum in water zone: build water first, derive oil from WOC
                auto wDensFn = [&](Scalar depth, Scalar press) {
                    return waterDensityFn(depth, press);
                };
                waterPressFunc = buildPressureTable(datumDepth, datumPressure,
                                                    wDensFn);
                waterBuilt = true;
                const Scalar oilPAtWoc = waterPressFunc.eval(wocDepth, true) + wocPc;
                oilPressFunc = buildPressureTable(wocDepth, oilPAtWoc, oilDensityFn);
                oilBuilt = true;
            } else {
                // Datum in oil zone (most common): build oil directly
                oilPressFunc = buildPressureTable(datumDepth, datumPressure,
                                                  oilDensityFn);
                oilBuilt = true;
            }
        }

        // Build gas pressure table if not yet built
        if (gas_active && !gasBuilt) {
            auto gasDensityFn = [&](Scalar depth, Scalar press) {
                return hydrocarbonDensity(depth, press, gasPhaseIdx);
            };
            if (oilBuilt) {
                const Scalar gasPAtGoc = oilPressFunc.eval(gocDepth, true) + gocPc;
                gasPressFunc = buildPressureTable(gocDepth, gasPAtGoc, gasDensityFn);
            } else {
                gasPressFunc = buildPressureTable(datumDepth, datumPressure,
                                                  gasDensityFn);
            }
            gasBuilt = true;
        }

        // Build water pressure table if not yet built
        if (water_active && !waterBuilt) {
            auto wDensFn = [&](Scalar depth, Scalar press) {
                return waterDensityFn(depth, press);
            };
            if (oilBuilt) {
                const Scalar wPAtWoc = oilPressFunc.eval(wocDepth, true) - wocPc;
                waterPressFunc = buildPressureTable(wocDepth, wPAtWoc, wDensFn);
            } else {
                waterPressFunc = buildPressureTable(datumDepth, datumPressure,
                                                    wDensFn);
            }
            waterBuilt = true;
        }

        // --- Populate initial fluid states ---
        zmf_initialization_ = true;

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& fs = initialFluidStates_[dofIdx];
            const Scalar depth = cellCenterDepths[dofIdx];

            // Temperature
            fs.setTemperature(tempVdFunc.eval(depth, true));

            // Phase pressures from the integrated tables
            if (oil_active)
                fs.setPressure(oilPhaseIdx, oilPressFunc.eval(depth, true));
            if (gas_active)
                fs.setPressure(gasPhaseIdx, gasPressFunc.eval(depth, true));
            if (water_active)
                fs.setPressure(waterPhaseIdx, waterPressFunc.eval(depth, true));

            // Saturations based on contact depths
            Scalar Sg = Scalar(0), So = Scalar(0), Sw = Scalar(0);
            if (gas_active && oil_active && depth <= gocDepth) {
                // Gas zone (above GOC)
                Sg = Scalar(1);
            } else if (water_active && oil_active && depth >= wocDepth) {
                // Water zone (below WOC)
                Sw = Scalar(1);
            } else {
                // Oil zone (between contacts) or the only active zone
                if (oil_active)
                    So = Scalar(1);
                else if (gas_active)
                    Sg = Scalar(1);
                else
                    Sw = Scalar(1);
            }

            if (gas_active) fs.setSaturation(gasPhaseIdx, Sg);
            if (oil_active) fs.setSaturation(oilPhaseIdx, So);
            if (water_active) fs.setSaturation(waterPhaseIdx, Sw);

            // Overall (total) mole fractions from ZMFVD or uniform default
            std::vector<Scalar> comp(numComponents);
            compositionAtDepth(depth, comp);
            for (unsigned c = 0; c < numComponents; ++c)
                fs.setMoleFraction(c, comp[c]);

            // Phase compositions (initial guess; the simulator will flash)
            for (unsigned c = 0; c < numComponents; ++c) {
                if (oil_active)
                    fs.setMoleFraction(oilPhaseIdx, c, comp[c]);
                if (gas_active)
                    fs.setMoleFraction(gasPhaseIdx, c, comp[c]);
            }
        }
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
