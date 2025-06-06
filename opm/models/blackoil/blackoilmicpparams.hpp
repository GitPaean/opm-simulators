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
/*!
 * \file
 *
 * \brief Contains the parameters required to extend the black-oil model by MICP.
 */
#ifndef OPM_BLACK_OIL_MICP_PARAMS_HPP
#define OPM_BLACK_OIL_MICP_PARAMS_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <vector>

namespace Opm {

#if HAVE_ECL_INPUT
class EclipseState;
#endif

//! \brief Struct holding the parameters for the BlackOilMICPModule class.
template<class Scalar>
struct BlackOilMICPParams
{
#if HAVE_ECL_INPUT
    template<bool enableMICP>
    void initFromState(const EclipseState& eclState);
#endif

    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    std::vector<TabulatedFunction> permfactTable_;

    std::vector<Scalar> densityBiofilm_;
    std::vector<Scalar> densityCalcite_;
    std::vector<Scalar> detachmentRate_;
    std::vector<Scalar> detachmentExponent_;
    std::vector<Scalar> halfVelocityOxygen_;
    std::vector<Scalar> halfVelocityUrea_;
    std::vector<Scalar> maximumGrowthRate_;
    std::vector<Scalar> maximumUreaUtilization_;
    std::vector<Scalar> microbialAttachmentRate_;
    std::vector<Scalar> microbialDeathRate_;
    std::vector<Scalar> oxygenConsumptionFactor_;
    std::vector<Scalar> yieldGrowthCoefficient_;
    std::vector<Scalar> yieldUreaToCalciteCoefficient_;
    std::vector<Scalar> microbialDiffusion_;
    std::vector<Scalar> oxygenDiffusion_;
    std::vector<Scalar> ureaDiffusion_;
};

} // namespace Opm

#endif // OPM_BLACK_OIL_MICP_PARAMS_HPP
