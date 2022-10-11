/*
  Copyright 2021 Equinor ASA.


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

#ifndef OPM_PERFDATA_HEADER_INCLUDED
#define OPM_PERFDATA_HEADER_INCLUDED

#include <vector>

namespace Opm
{

class PerfData
{
private:
    bool injector;

public:
    PerfData(std::size_t num_perf, double pressure_first_connection_, bool injector_, std::size_t num_phases);
    std::size_t size() const;
    bool empty() const;
    bool try_assign(const PerfData& other);


    double pressure_first_connection;
    std::vector<double> pressure;
    std::vector<double> rates;
    std::vector<double> phase_rates;
    std::vector<double> solvent_rates;
    std::vector<double> polymer_rates;
    std::vector<double> brine_rates;
    std::vector<double> prod_index;
    std::vector<double> micp_rates;

    std::vector<std::size_t> cell_index;
    std::vector<double> connection_transmissibility_factor;
    // for keyword WINJMULT, only useful for injectors with CIRR mode
    std::vector<double> inj_multipler;
    std::vector<int> satnum_id;
    std::vector<std::size_t> ecl_index;


    // The water_throughput, skin_pressure and water_velocity variables are only
    // used for injectors to check the injectivity.
    std::vector<double> water_throughput;
    std::vector<double> skin_pressure;
    std::vector<double> water_velocity;
};

} // namespace Opm

#endif // OPM_PERFORATIONDATA_HEADER_INCLUDED
