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

#include <algorithm>
#include <functional>
#include <vector>

namespace Opm::Parameters {

// Reproduce the reference simulator's normalization of an explicitly
// initialized water saturation in compositional runs (see
// FlowProblemComp::readExplicitInitialCondition_). Off by default: the plain
// reading of SWAT is used unless this compatibility behaviour is requested.
struct CompatExplicitSwatInit { static constexpr bool value = false; };

} // namespace Opm::Parameters

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

        Parameters::Register<Parameters::CompatExplicitSwatInit>
            ("Normalize an explicitly initialized water saturation the way the "
             "reference simulator does: hydrocarbon moles are seeded from the "
             "unflashed overall composition filling (1 - SWAT) of the pore "
             "space, and after the flash all phase volumes are rescaled to "
             "fill the pore volume, so the water saturation deviates from the "
             "input SWAT by the flash volume change.");
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

            this->eclWriter_->writeOutput(std::move(localCellData), isSubStep,
                                          this->simulator().vanguard().schedule()
                                         .exitStatus().has_value());
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

    void setSubStepReport(const SimulatorReportSingle& report)
    { return eclWriter_->setSubStepReport(report); }

    void setSimulationReport(const SimulatorReport& report)
    { return eclWriter_->setSimulationReport(report); }

    void finalizeOutput()
    {
        OPM_TIMEBLOCK(finalizeOutput);
        eclWriter_.reset();
    }

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
        throw std::logic_error("Equilibration is not supported by compositional modeling yet");
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

        if (water_active && zmf_initialization_ &&
            Parameters::Get<Parameters::CompatExplicitSwatInit>())
        {
            compatExplicitSwatInit_();
        }
    }

    // The reference simulator does not use an explicitly initialized SWAT as
    // the water saturation directly. It seeds the cell's hydrocarbon amount
    // as if (1 - SWAT) of the pore space held the *unflashed* overall
    // composition at its single-phase EOS molar volume, flashes, and then
    // rescales all phase volumes by a common factor to fill the pore volume.
    // The water saturation therefore deviates from the input by the flash
    // volume change:
    //
    //   Sw / (1 - Sw) = [SWAT / (1 - SWAT)] * Vm_single(z) / Vm_flash(z)
    //
    // This was reverse engineered from the reference outputs (matched to
    // ~1e-5 in Sw on fluids with and without binary interaction
    // coefficients) and is reproduced here for comparison runs only.
    void compatExplicitSwatInit_()
    {
        using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;
        const auto eos_type = getEosType();
        constexpr Scalar flash_tolerance = 1e-8;

        for (auto& fs : initialFluidStates_) {
            const Scalar sw_in = fs.saturation(FluidSystem::waterPhaseIdx);
            if (sw_in <= 0.0 || sw_in >= 1.0) {
                continue;
            }

            InitialFluidState flash_fs;
            flash_fs.setTemperature(fs.temperature(0));
            flash_fs.setPressure(FluidSystem::oilPhaseIdx, fs.pressure(FluidSystem::oilPhaseIdx));
            flash_fs.setPressure(FluidSystem::gasPhaseIdx, fs.pressure(FluidSystem::gasPhaseIdx));
            for (unsigned c = 0; c < numComponents; ++c) {
                flash_fs.setMoleFraction(c, fs.moleFraction(c));
                flash_fs.setKvalue(c, flash_fs.wilsonK_(c));
            }
            flash_fs.setLvalue(-1.0);
            const bool single_phase =
                FlashSolver::flash_solve_scalar_(flash_fs, "ssi", flash_tolerance, eos_type);
            if (single_phase) {
                continue;   // no volume change, SWAT stands as given
            }

            typename FluidSystem::template ParameterCache<Scalar> param_cache(eos_type);

            // molar volume of the unflashed overall composition (the EOS root
            // the liquid-phase selection yields; unique on the states tested)
            InitialFluidState single_fs(flash_fs);
            for (unsigned c = 0; c < numComponents; ++c) {
                single_fs.setMoleFraction(FluidSystem::oilPhaseIdx, c, fs.moleFraction(c));
            }
            param_cache.updatePhase(single_fs, FluidSystem::oilPhaseIdx);
            const Scalar vm_single = param_cache.molarVolume(FluidSystem::oilPhaseIdx);

            param_cache.updatePhase(flash_fs, FluidSystem::oilPhaseIdx);
            const Scalar vm_oil = param_cache.molarVolume(FluidSystem::oilPhaseIdx);
            param_cache.updatePhase(flash_fs, FluidSystem::gasPhaseIdx);
            const Scalar vm_gas = param_cache.molarVolume(FluidSystem::gasPhaseIdx);
            const Scalar L = flash_fs.L();
            const Scalar vm_flash = L * vm_oil + (1.0 - L) * vm_gas;

            const Scalar ratio = sw_in / (1.0 - sw_in) * vm_single / vm_flash;
            const Scalar sw = ratio / (1.0 + ratio);
            const Scalar hc_scale = (1.0 - sw) / (1.0 - sw_in);

            fs.setSaturation(FluidSystem::waterPhaseIdx, sw);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                fs.setSaturation(FluidSystem::gasPhaseIdx,
                                 fs.saturation(FluidSystem::gasPhaseIdx) * hc_scale);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                fs.setSaturation(FluidSystem::oilPhaseIdx,
                                 fs.saturation(FluidSystem::oilPhaseIdx) * hc_scale);
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
