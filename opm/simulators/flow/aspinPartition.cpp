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
#include "config.h"

#include <opm/simulators/flow/aspinPartition.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>

// TODO: cleaning the header files from the opm-grid
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>

namespace Opm
{


namespace
{

    // Read from file, containing one number per cell, from [0, ... , num_domains - 1].
    std::pair<std::vector<int>, int> partitionCellsFromFile([[maybe_unused]] const int num_cells)
    {
        // TODO: refactor to make more flexible.
        // Read file into single vector.
        const std::string filename = "partition.txt";
        std::ifstream is(filename);
        const std::vector<int> cellpart{std::istream_iterator<int>(is), std::istream_iterator<int>()};
        if (cellpart.size() != size_t(num_cells)) {
            auto msg = fmt::format("Partition file contains {} entries, but there are {} cells.",
                                   cellpart.size(), num_cells);
            throw std::runtime_error(msg);
        }

        // Create and return the output domain vector.
        const int num_domains = (*std::max_element(cellpart.begin(), cellpart.end())) + 1;
        return { cellpart, num_domains };
    }


    // Trivially simple partitioner
    std::pair<std::vector<int>, int> partitionCellsSimple(const int num_cells, const int num_domains)
    {
        // Build the partitions.
        std::vector<int> bounds(num_domains + 1);
        bounds[0] = 0;
        const int dom_sz = num_cells / num_domains;
        for (int i = 1; i < num_domains; ++i) {
            bounds[i] = dom_sz * i;
        }
        bounds[num_domains] = num_cells;
        std::vector<int> part(num_cells);
        for (int i = 0; i < num_domains; ++i) {
            std::fill(part.begin() + bounds[i], part.begin() + bounds[i + 1], i);
        }
        return { part, num_domains };
    }

    std::pair<std::vector<int>, int> partitionWithZoltan(const int num_cells, const int num_domains,
                                                         const Dune::CpGrid& grid, const std::vector<Well>& wells) {
//        return std::pair<std::vector<int>, int>{};
        auto cc = Dune::MPIHelper::getCollectiveCommunication();
        struct Zoltan_Struct *zz;
        zz = Zoltan_Create(MPI_COMM_WORLD);
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, grid);


        // using 0 uniformEdgeWgt for now
        std::vector<int> computedCellPart;
        std::vector<std::pair<std::string,bool>> wells_on_proc;
        std::vector<std::tuple<int,int,char>> exportList;
        std::vector<std::tuple<int,int,char,int>> importList;
        std::tie(computedCellPart, wells_on_proc, exportList, importList) =
                Dune::cpgrid::zoltanSerialGraphPartitionGridOnRoot(grid, &wells, nullptr, cc, Dune::uniformEdgeWgt, 0, 1.1, false);

        // return std::pair<std::vector<int>, int>{};
        return std::make_pair(computedCellPart, 2);


        // not handling transmissibilites for now, assuming homogenous transmissibilities
        // TODO: with this way, looking for the needed partition information, which process the grid cell is allocated at.
        // /// @return A tuple consisting of a vector that contains for each local cell of the original grid the
        /////         the number of the process that owns it after repartitioning,
        /////         a set of names of wells that should be defunct in a parallel
        /////         simulation, vector containing information for each exported cell (global id
        /////         of cell, process id to send to, attribute there), and a vector containing
        /////         information for each imported cell (global index, process id that sends, attribute here, local index
        /////         here)

        /* std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
                std::vector<std::tuple<int,int,char> >,
                std::vector<std::tuple<int,int,char,int> >  >
        zoltanSerialGraphPartitionGridOnRoot(const CpGrid& grid,
                                             const std::vector<OpmWellType> * wells,
                                             const double* transmissibilities,
                                             const CollectiveCommunication<MPI_Comm>& cc,
                                             EdgeWeightMethod edgeWeightsMethod, int root,
                                             const double zoltanImbalanceTol,
                                             bool allowDistributedWells); */
        // Basically, we need the following three element for here,
        // CpGrid, Wells, and transmissibilities
        /* for CpGrid,
         * for wells, schedule.getWellsatEnd(), it needs all the wells
         *
         * for transmissibilities
         * faceTrans.resize(numFaces, 0.0);
            ElementMapper elemMapper(gridv, Dune::mcmgElementLayout()); */
            // auto elemIt = gridView.template begin</*codim=*/0>();
        // const auto& elemEndIt = gridView.template end</*codim=*/0>();
        /* for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;
            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                const auto& is = *isIt;
                if (!is.neighbor())
                    continue;

                unsigned I = elemMapper.index(is.inside());
                unsigned J = elemMapper.index(is.outside());

                // FIXME (?): this is not portable!
                unsigned faceIdx = is.id();

                faceTrans[faceIdx] = this->getTransmissibility(I,J);
         *
         *
         */
    }

} // anonymous namespace


std::pair<std::vector<int>, int>
partitionCells(const int num_cells, const Dune::CpGrid& grid, const std::vector<Well>& wells)
{
    // const std::string method = "simple";
    const std::string method = "zoltan";
    if (method == "simple") {
        return partitionCellsSimple(num_cells, 25);
    } else if (method == "file") {
        return partitionCellsFromFile(num_cells);
    } else if (method == "zoltan") {
        return partitionWithZoltan(num_cells, 25, grid, wells);
    } else {
        return {};
    }
}

} // namespace Opm
