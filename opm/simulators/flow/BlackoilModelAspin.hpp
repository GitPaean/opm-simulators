/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS
  Copyright 2021 SINTEF Digital

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

#ifndef OPM_BLACKOILMODELASPIN_HEADER_INCLUDED
#define OPM_BLACKOILMODELASPIN_HEADER_INCLUDED


#include <opm/simulators/flow/BlackoilModelEbos.hpp>

namespace Opm {
    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    /// This variant implements an Additive Schwarz Preconditioned Inexact Newton (ASPIN) nonlinear solution scheme.
    template <class TypeTag>
    class BlackoilModelAspin
    {
    public:
        using BaseModel = BlackoilModelEbos<TypeTag>;
        using Simulator       = typename BaseModel::Simulator;
        using ModelParameters = typename BaseModel::ModelParameters;
        using WellState       = typename BaseModel::WellState;

        // ---------  Public methods  ---------

        /// Construct the model.
        BlackoilModelAspin(Simulator& ebosSimulator,
                           const ModelParameters& param,
                           BlackoilWellModel<TypeTag>& well_model,
                           const bool terminal_output)
            : base_model_(ebosSimulator, param, well_model, terminal_output)
        {
        }

        void beginReportStep()
        {
            base_model_.beginReportStep();
        }

        void endReportStep()
        {
            base_model_.endReportStep();
        }

        const auto& wellModel() const
        {
            return base_model_.wellModel();
        }

        auto& wellModel()
        {
            return base_model_.wellModel();
        }

        auto numPhases() const
        {
            return base_model_.numPhases();
        }

        Simulator& ebosSimulator()
        {
            return base_model_.ebosSimulator();
        }

        const auto& failureReport() const
        {
            return base_model_.failureReport();
        }

        const auto& stepReports() const
        {
            return base_model_.stepReports();
        }

        auto computeFluidInPlace(const std::vector<int>& fipnum) const
        {
            // TODO: this is weird, check it.
            return base_model_.computeFluidInPlace(fipnum);
        }

        double relativeChange() const
        {
            return base_model_.relativeChange();
        }

        SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer)
        {
            return base_model_.prepareStep(timer);
        }

        SimulatorReportSingle afterStep(const SimulatorTimerInterface& timer)
        {
            return base_model_.afterStep(timer);
        }

        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const int iteration,
                                                 const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
        {
            SimulatorReportSingle report;
            failureReport_ = SimulatorReportSingle();
            Dune::Timer perfTimer;

            perfTimer.start();
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
                convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
                convergence_reports_.back().report.reserve(11);
            }

            report.total_linearizations = 1;

            try {
                report += base_model_.assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }

            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration,residual_norms);
                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(Opm::NumericalIssue, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW(Opm::NumericalIssue, "Too large residual found!");
                }
            }
            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                const int nc = UgGridHelpers::numCells(grid_);
                BVector x(nc);

                // apply the Schur compliment of the well model to the reservoir linearized
                // equations
                wellModel().linearize(ebosSimulator().model().linearizer().jacobian(),
                                      ebosSimulator().model().linearizer().residual());

                // Solve the linear system.
                linear_solve_setup_time_ = 0.0;
                try {
                    solveJacobianSystem(x);
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();
                }
                catch (...) {
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();

                    failureReport_ += report;
                    throw; // re-throw up
                }

                perfTimer.reset();
                perfTimer.start();

                // handling well state update before oscillation treatment is a decision based
                // on observation to avoid some big performance degeneration under some circumstances.
                // there is no theorectical explanation which way is better for sure.
                wellModel().postSolve(x);

                if (param_.use_update_stabilization_) {
                    // Stabilize the nonlinear update.
                    bool isOscillate = false;
                    bool isStagnate = false;
                    nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                    if (isOscillate) {
                        current_relaxation_ -= nonlinear_solver.relaxIncrement();
                        current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                        if (terminalOutputEnabled()) {
                            std::string msg = "    Oscillating behavior detected: Relaxation set to "
                                    + std::to_string(current_relaxation_);
                            OpmLog::info(msg);
                        }
                    }
                    nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
                }

                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateSolution(x);

                report.update_time += perfTimer.stop();
            }

            return report;
        }

    private:
        BaseModel base_model_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELASPIN_HEADER_INCLUDED
