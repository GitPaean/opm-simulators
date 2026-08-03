/*
  Copyright 2026, SINTEF Digital

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
#ifndef OPM_COMPOSITIONAL_CASE_CHECK_HPP
#define OPM_COMPOSITIONAL_CASE_CHECK_HPP

#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <fmt/format.h>

#include <cstdio>
#include <string>

namespace Opm
{

/*!
 * \brief Report a deck the compositional simulators cannot run.
 *
 * CO2STORE and H2STORE are black-oil options, and Runspec::compositional()
 * excludes them even when COMPS is present, while a deck without COMPS has no
 * components to solve for.  None of these end up with a usable compositional
 * fluid system, so say so while the deck is the only thing that has been read
 * rather than failing deep inside problem construction.
 *
 * \param[in] runspec RUNSPEC section of the deck.
 *
 * \param[in] deckFile Deck the run was started with, for the diagnostic.
 *
 * \return Whether or not the run may proceed.
 */
inline bool
isCompositionalCase(const Runspec& runspec, const std::string& deckFile)
{
    if (runspec.compositional()) {
        return true;
    }

    // Runspec::compositional() is (numComps() > 0) && !co2Storage() &&
    // !h2Storage(), so these three cases are exhaustive here.
    const auto* reason = runspec.co2Storage() ? "it is a CO2STORE run"
        : runspec.h2Storage()                 ? "it is an H2STORE run"
                                              : "it has no COMPS keyword";

    fmt::print(stderr,
               "{} is not a compositional case: {}.\n"
               "This simulator runs compositional cases only - "
               "use flow for this deck.\n",
               deckFile,
               reason);

    return false;
}

} // namespace Opm

#endif // OPM_COMPOSITIONAL_CASE_CHECK_HPP
