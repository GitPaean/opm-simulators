/*
  Copyright 2021 Total SE

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

#include <opm/simulators/flow/aspinPartition.hpp>
#include <numeric>

namespace Opm
{

std::vector<std::vector<int>> partitionCells(const int num_cells)
{
    // TODO: replace with non-trivial version...
    const int num_domains = 4;

    // Local lambda to return cell indices of a domain.
    auto createCells = [num_cells](int dom_index) {
        std::vector<int> cells;
        const int dom_sz = num_cells/num_domains;
        const int start = dom_index * dom_sz;
        if (dom_index < num_domains - 1) {
            // Not the last domain
            cells.resize(dom_sz);
        } else {
            // The last domain
            cells.resize(num_cells - start);
        }
        std::iota(cells.begin(), cells.end(), start);
        return cells;
    };

    // Build the partitions.
    std::vector<std::vector<int>> part(num_domains);
    for (int i = 0; i < num_domains; ++i) {
        part[i] = createCells(i);
    }
    return part;
}

} // namespace Opm


