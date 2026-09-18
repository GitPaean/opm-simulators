/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_WELL_RATE_ALLOCATION_HPP
#define OPM_WELL_RATE_ALLOCATION_HPP

#include <cmath>
#include <optional>

namespace Opm {

template<class Scalar>
struct ReportedPhaseSplit
{
    Scalar free;
    Scalar solution;
};

// Signed rates: production is negative. rateTolerance has rate units.
// This is only a reporting allocation for production without crossflow.
template<class Scalar>
std::optional<ReportedPhaseSplit<Scalar>>
allocateProductionRate(Scalar total, Scalar free, Scalar solution,
                       Scalar rateTolerance, bool hasCrossflow)
{
    if (!std::isfinite(total) || !std::isfinite(free) ||
        !std::isfinite(solution) || !std::isfinite(rateTolerance) ||
        rateTolerance < Scalar{0}) {
        return std::nullopt;
    }

    if (total == Scalar{0}) {
        return ReportedPhaseSplit<Scalar>{Scalar{0}, Scalar{0}};
    }

    if (hasCrossflow || total > Scalar{0} ||
        free > Scalar{0} || solution > Scalar{0}) {
        return std::nullopt;
    }

    const Scalar sum = free + solution;
    if (!std::isfinite(sum) || -sum <= rateTolerance) {
        return std::nullopt;
    }

    const Scalar reportedFree = total * (free / sum);
    return ReportedPhaseSplit<Scalar>{reportedFree, total - reportedFree};
}

} // namespace Opm

#endif // OPM_WELL_RATE_ALLOCATION_HPP
