/*
  Copyright 2023 Equinor ASA.

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

#include <config.h>

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <iostream>
#include <cstdlib>

Opm::ParallelPAvgDynamicSourceData::
ParallelPAvgDynamicSourceData(const Parallel::Communication&  comm,
                              const std::vector<std::size_t>& sourceLocations,
                              GlobalToLocal                   localCellIdx)
    : PAvgDynamicSourceData { sourceLocations }
    , comm_                 { comm }
{
    std::cout << " ParallelPAvgDynamicSourceData sourlocatoins.size() " << sourceLocations.size() << std::endl;
    this->finaliseConstruction(sourceLocations, std::move(localCellIdx));
}

void Opm::ParallelPAvgDynamicSourceData::setToZero()
{
    std::fill_n(this->localSrc_.begin(), this->localSrc_.size(), 0.0);
}

void
Opm::ParallelPAvgDynamicSourceData::
reconstruct(const std::vector<std::size_t>& sourceLocations,
            GlobalToLocal                   localCellIdx)
{
    PAvgDynamicSourceData::reconstruct(sourceLocations); // Reconstruct base
    this->finaliseConstruction(sourceLocations, std::move(localCellIdx));
}

void Opm::ParallelPAvgDynamicSourceData::collectLocalSources(Evaluator eval)
{
    auto localIx = std::size_t{0};

    for (const auto& location : this->locations_) {
        eval(location.cell, this->localSourceTerm(localIx++));
    }
}

void Opm::ParallelPAvgDynamicSourceData::synchroniseSources()
{
    this->comm_.get()
        .allgatherv(this->localSrc_.data(),       // Input (from)
                    static_cast<int>(this->localSrc_.size()),
                    this->src_.data(),            // Output (to)
                    this->allSizes_.data(),       // #elements per rank
                    this->startPointers_.data()); // Output offsets
}

std::vector<double>::size_type
Opm::ParallelPAvgDynamicSourceData::
storageIndex(const std::vector<double>::size_type elemIndex) const
{
    return this->storageIndex_[elemIndex];
}

void
Opm::ParallelPAvgDynamicSourceData::
finaliseConstruction(const std::vector<std::size_t>& sourceLocations,
                     GlobalToLocal                   localCellIdx)
{
    std::cout << " finaliseConstruction sourceLocations.size() " << sourceLocations.size() << std::endl;
    auto ix = std::size_t{0};

    this->locations_.clear();

    for (const auto& location : sourceLocations) {
        if (const auto cell = localCellIdx(location); cell >= 0) {
            this->locations_.push_back({ ix, cell });
        }

        ix += 1;
    }
    const bool output_50 = this->locations_.size() == 50;
    if (output_50) {
        std::cout << " outputting the locations_ construction process " << std::endl;
        for (const auto& location : sourceLocations) {
            const auto cell = localCellIdx(location);
            std::cout << " location " << location << " cell " << cell << std::endl;
        }
        size_t count = 0;
        std::cout << " outputting the one with cell > 0" << std::endl;
        for (const auto& location : sourceLocations) {
            if (const auto cell = localCellIdx(location); cell >= 0) {
                std::cout << " location " << location << " cell " << cell << std::endl;
            }
            count++;
        }
        std::cout << " there are " << count << " cell >=0 " << std::endl;
    }

    this->localSrc_.assign(numSpanItems() * this->locations_.size(), 0.0);

    this->defineCommunication();
}

Opm::PAvgDynamicSourceData::SourceDataSpan<double>
Opm::ParallelPAvgDynamicSourceData::localSourceTerm(const std::size_t localIx)
{
    return this->sourceTerm(localIx, this->localSrc_);
}

void Opm::ParallelPAvgDynamicSourceData::defineCommunication()
{
    std::cout << "defineCommunication" << std::endl;
    // 1) Determine origins/owning ranks for all source terms.
    const bool output_50 = this->locations_.size() == 50;
    if (output_50) {
        std::cout << " output this->locations " << std::endl;
        for (size_t i = 0; i < this->locations_.size(); ++i) {
            std::cout << " i " << i << " ix " << this->locations_[i].ix << " cell " << this->locations_[i].cell << std::endl;
        }
    }
    auto ixVec = std::vector<std::size_t>(this->locations_.size());
    std::transform(this->locations_.begin(), this->locations_.end(),
                   ixVec.begin(),
                   [](const auto& location) { return location.ix; });

    constexpr auto numItems = numSpanItems();
    if (output_50) {
        std::cout << " output ixVec " << std::endl;
        for (size_t i = 0; i < ixVec.size(); ++i) {
            std::cout << " " << ixVec[i];
            if ((i+1) % 10 == 0) {
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }

    const auto& [allIndices, allIxStart] = allGatherv(ixVec, this->comm_.get());
    if (output_50) {
        std::cout << " output allIndices " << std::endl;
        for (size_t i = 0; i < allIndices.size(); ++i) {
            std::cout << " " << allIndices[i];
            if ((i + 1) % 10 == 0) {
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }

    // -----------------------------------------------------------------------

    // 2) Determine starting pointers/offsets/displacements for received
    //    basic elements from each rank.  There are 'numItems' basic data
    //    elements for each source term.
    this->startPointers_.resize(allIxStart.size());
    std::transform(allIxStart.begin(), allIxStart.end(),
                   this->startPointers_.begin(),
                   [](const int start)
                   {
                       return numItems * start;
                   });

    // -----------------------------------------------------------------------

    // 3) Determine number of basic data elements to receive from each rank.
    this->allSizes_.resize(allIxStart.size() - 1);
    std::adjacent_difference(this->startPointers_.begin() + 1,
                             this->startPointers_.end(),
                             this->allSizes_.begin());

    // -----------------------------------------------------------------------

    // 4) Build translation mapping from source term element indices to
    //    storage indices.
    this->storageIndex_.resize(allIndices.size());
    bool potential_invalid_write = false;
    auto storageIx = std::vector<double>::size_type{0};
    for (const auto& elemIndex : allIndices) {
        if (elemIndex >= this->storageIndex_.size()) {
            std::cout << " elemIndex is " << elemIndex << " storageIndex_.size() is " << this->storageIndex_.size() << " storageIx " << storageIx << std::endl;
            potential_invalid_write = true;
        }
        this->storageIndex_[elemIndex] = storageIx++;
    }

    if (potential_invalid_write) {
        std::cout << " potential_invalid_write happens, outputting allIndices now " << std::endl;
        for (size_t i = 0; i < allIndices.size(); ++i) {
            std::cout << " " << allIndices[i];
            if ( (i + 1)%5 == 0) {
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
        std::cout << " output this->locations_ " << std::endl;
        for (size_t i = 0; i < this->locations_.size(); ++i) {
            std::cout << " " << this->locations_[i].ix << " " << this->locations_[i].cell << std::endl;
        }
        exit(100);
    }

}
