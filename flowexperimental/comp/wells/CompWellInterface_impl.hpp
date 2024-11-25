/*
  Copyright 2024, SINTEF Digital

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

#include "CompConnectionData.hpp"

#include <string>

namespace Opm
{
template <typename TypeTag>
CompWellInterface<TypeTag>::
CompWellInterface(const Well& well,
                  const int index_of_well,
                  const std::vector<CompConnectionData<Scalar> >& well_connection_data)
    : well_ecl_(well)
    , index_of_well_(index_of_well)
    , reference_depth_(well.getRefDepth())
{
    number_of_connection_ = well_connection_data.size();
    {
        well_cells_.resize(number_of_connection_);
        well_index_.resize(number_of_connection_);
        saturation_table_number_.resize(number_of_connection_, 0);
        int connection_idx = 0;
        for (const auto& connection_data : well_connection_data) {
            well_cells_[connection_idx] = connection_data.cell_index;
            well_index_[connection_idx] = connection_data.connection_transmissibility_factor;
            saturation_table_number_[connection_idx] = connection_data.satnum_id;
            ++connection_idx;
        }
        // TODO: saturation_table_number
    }

}

template <typename TypeTag>
void
CompWellInterface<TypeTag>::
init() {
    // more things to add here
}

} // end of namespace Opm