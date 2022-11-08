/*
  Copyright Equinor ASA 2021

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
#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <algorithm>
#include <stdexcept>

#include <opm/simulators/wells/SegmentState.hpp>

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <sstream>

namespace Opm
{

namespace {
std::vector<int> make_segment_number( const WellSegments& segments ) {
    std::vector<int> segment_number;
    std::transform(segments.begin(), segments.end(), std::back_insert_iterator(segment_number), [](const Segment& segment) { return segment.segmentNumber(); });
    return segment_number;
}

}

SegmentState::SegmentState(int num_phases, const WellSegments& segments) :
    rates(segments.size() * num_phases),
    pressure(segments.size()),
    pressure_drop_friction(segments.size()),
    pressure_drop_hydrostatic(segments.size()),
    pressure_drop_accel(segments.size()),
    m_segment_number(make_segment_number(segments))
{
}

double SegmentState::pressure_drop(std::size_t index) const {
    return this->pressure_drop_friction[index] + this->pressure_drop_hydrostatic[index] + this->pressure_drop_accel[index];
}


bool SegmentState::empty() const {
    return this->rates.empty();
}

std::size_t SegmentState::size() const {
    return this->pressure.size();
}


void SegmentState::scale_pressure(const double bhp) {
    if (this->empty())
        throw std::logic_error("Tried to pressure scale empty SegmentState");

    const auto pressure_change = bhp - this->pressure[0];

    std::transform(this->pressure.begin(),
                   this->pressure.end(),
                   this->pressure.begin(),
                   [pressure_change] (const double& p) { return p + pressure_change;});
}

const std::vector<int>& SegmentState::segment_number() const {
    return this->m_segment_number;
}

std::string SegmentState::output() const {
    std::stringstream ss;
    const size_t number_phases = this->rates.size() / this->pressure.size();
    const size_t nseg = this->size();
    ss << " outputting the segment information, there are " << nseg << " segments, and " << number_phases << " phases, rates and pressures are :" << std::endl;
    for (size_t seg = 0; seg < nseg; ++seg) {
        ss  << " seg " << seg << " rates ";
        for (size_t p = 0; p < number_phases; ++p) {
            ss << " " << this->rates[number_phases * seg + p] * 86400.;
        }
        ss << " pressure " << this->pressure[seg] / 1.e5
           << " pressure_drop_friction " << this->pressure_drop_friction[seg] / 1.e5 << " pressure_drop_hydrostatic " << this->pressure_drop_hydrostatic[seg]/1.e5
           << " pressure_drop_accel " << this->pressure_drop_accel[seg] / 1.e5 << std::endl;
    }
    return ss.str();
}

}
