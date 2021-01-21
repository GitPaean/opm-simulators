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
            return base_model_.nonlinearIteration(iteration, timer, nonlinear_solver);
        }

    private:
        BaseModel base_model_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELASPIN_HEADER_INCLUDED
