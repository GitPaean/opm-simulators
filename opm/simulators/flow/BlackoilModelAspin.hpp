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

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModelAspin(Simulator& ebosSimulator,
                           const ModelParameters& param,
                           BlackoilWellModel<TypeTag>& well_model,
                           const bool terminal_output)
            : base_model_(ebosSimulator, param, well_model, terminal_output)
        {
        }
    private:
        BaseModel base_model_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELASPIN_HEADER_INCLUDED
