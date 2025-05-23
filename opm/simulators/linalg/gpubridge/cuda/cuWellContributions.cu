/*
  Copyright 2020 Equinor ASA

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


#include <config.h> // CMake
#include <opm/simulators/linalg/gpubridge/cuda/cuWellContributions.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/linalg/gpubridge/MultisegmentWellContribution.hpp>
#include <opm/simulators/linalg/gpubridge/cuda/cuda_header.hpp>

#include <cuda_runtime.h>

namespace Opm
{

// apply WellContributions using y -= C^T * (D^-1 * (B * x))
template<class Scalar>
__global__ void apply_well_contributions(const Scalar* __restrict__ Cnnzs,
                                         const Scalar* __restrict__ Dnnzs,
                                         const Scalar* __restrict__ Bnnzs,
                                         const int* __restrict__ Ccols,
                                         const int* __restrict__ Bcols,
                                         const Scalar* __restrict__ x,
                                         Scalar* __restrict__ y,
                                         const int dim,
                                         const int dim_wells,
                                         const unsigned int * __restrict__ val_pointers)
{
    const int idx_b = blockIdx.x;
    const int idx_t = threadIdx.x;
    const unsigned int val_size = val_pointers[idx_b + 1] - val_pointers[idx_b];

    const int vals_per_block = dim * dim_wells;        // 12
    const int num_active_threads = (32 / vals_per_block) * vals_per_block; // 24
    const int num_blocks_per_warp = 32 / vals_per_block; // 2
    const int lane = idx_t % 32;
    const int c = lane % dim;                           // col in block
    const int r = (lane / dim) % dim_wells;             // row in block

    extern __shared__ unsigned char smem[];
    Scalar* __restrict__ z1 = reinterpret_cast<Scalar*>(smem);
    Scalar* __restrict__ z2 = z1 + dim_wells;

    if (idx_t < dim_wells) {
        z1[idx_t] = 0.0;
    }

    __syncthreads();

    // z1 = B * x
    if (idx_t < num_active_threads) {
        // multiply all blocks with x
        Scalar temp = 0.0;
        int b = idx_t / vals_per_block + val_pointers[idx_b];       // block id, val_size indicates number of blocks
        while (b < val_size + val_pointers[idx_b]) {
            int colIdx = Bcols[b];
            temp += Bnnzs[b * dim * dim_wells + r * dim + c] * x[colIdx * dim + c];
            b += num_blocks_per_warp;
        }

        // merge all blocks into 1 dim*dim_wells block
        // since 3*4 blocks has give 2 parallel blocks, do not use a loop
        // 0x00ffffff contains 24 ones, representing the two blocks that are added
        // block 1:     block 2:
        //  0  1  2     12 13 14
        //  3  4  5     15 16 17
        //  6  7  8     18 19 20
        //  9 10 11     21 22 23
        // thread i will hold the sum of thread i and i + vals_per_block
        temp += __shfl_down_sync(0x00ffffff, temp, dim * dim_wells);

        // merge all (dim) columns of 1 block, results in a single 1*dim_wells vector, which is used to multiply with invD
        if (idx_t < vals_per_block) {
            // should be a loop as well, now only works for dim == 3
            if (c == 0 || c == 2) {temp += __shfl_down_sync(0x00000B6D, temp, 2);} // add col 2 to col 0
            if (c == 0 || c == 1) {temp += __shfl_down_sync(0x000006DB, temp, 1);} // add col 1 to col 0
        }

        // write 1*dim_wells vector to gmem, could be replaced with shfl broadcast to remove z1 altogether
        if (c == 0 && idx_t < vals_per_block) {
            z1[r] = temp;
        }
    }

    __syncthreads();

    // z2 = D^-1 * B * x = D^-1 * z1
    if (idx_t < dim_wells) {
        Scalar temp = 0.0;
        for (int c = 0; c < dim_wells; ++c) {
            temp += Dnnzs[idx_b * dim_wells * dim_wells + idx_t * dim_wells + c] * z1[c];
        }
        z2[idx_t] = temp;
    }

    __syncthreads();

    // y -= C^T * D^-1 * B * x
    // use dim * val_size threads, each block is assigned 'dim' threads
    if (idx_t < dim * val_size) {
        Scalar temp = 0.0;
        int b = idx_t / dim + val_pointers[idx_b];
        int cc = idx_t % dim;
        int colIdx = Ccols[b];
        for (unsigned int c = 0; c < dim_wells; ++c) {
            temp += Cnnzs[b * dim * dim_wells + c * dim + cc] * z2[c];
        }
        y[colIdx * dim + cc] -= temp;
    }
}

template<class Scalar>
WellContributionsCuda<Scalar>::~WellContributionsCuda()
{
    // delete data for StandardWell
    if (this->num_std_wells > 0) {
        cudaFree(d_Cnnzs);
        cudaFree(d_Dnnzs);
        cudaFree(d_Bnnzs);
        cudaFree(d_Ccols);
        cudaFree(d_Bcols);
        cudaFree(d_val_pointers);
    }

    if (this->num_ms_wells > 0 && h_x) {
        cudaFreeHost(h_x);
        cudaFreeHost(h_y);
        h_x = h_y = nullptr; // Mark as free for constructor
    }
}

template<class Scalar>
void WellContributionsCuda<Scalar>::APIalloc()
{
    cudaMalloc((void**)&d_Cnnzs,
               sizeof(Scalar) * this->num_blocks * this->dim * this->dim_wells);
    cudaMalloc((void**)&d_Dnnzs,
               sizeof(Scalar) * this->num_std_wells * this->dim_wells * this->dim_wells);
    cudaMalloc((void**)&d_Bnnzs,
               sizeof(Scalar) * this->num_blocks * this->dim * this->dim_wells);
    cudaMalloc((void**)&d_Ccols, sizeof(int) * this->num_blocks);
    cudaMalloc((void**)&d_Bcols, sizeof(int) * this->num_blocks);
    cudaMalloc((void**)&this->d_val_pointers, sizeof(unsigned int) * (this->num_std_wells + 1));
    cudaCheckLastError("apply_gpu malloc failed");
}

// Apply the WellContributions, similar to StandardWell::apply()
// y -= (C^T *(D^-1*(   B*x)))
template<class Scalar>
void WellContributionsCuda<Scalar>::apply(Scalar* d_x, Scalar* d_y)
{
    // apply MultisegmentWells

    // make sure the stream is empty if timing measurements are done
    cudaStreamSynchronize(stream);

    if (this->num_ms_wells > 0) {
        // allocate pinned memory on host if not yet done
        if (h_x == nullptr) {
            cudaMallocHost(&h_x, sizeof(Scalar) * this->N);
            cudaMallocHost(&h_y, sizeof(Scalar) * this->N);
        }

        // copy vectors x and y from GPU to CPU
        cudaMemcpyAsync(h_x, d_x, sizeof(Scalar) * this->N,
                        cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(h_y, d_y, sizeof(Scalar) * this->N,
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        // actually apply MultisegmentWells
        for (auto& well : this->multisegments) {
            well->apply(h_x, h_y);
        }

        // copy vector y from CPU to GPU
        cudaMemcpyAsync(d_y, h_y, sizeof(Scalar) * this->N,
                        cudaMemcpyHostToDevice, stream);
        cudaStreamSynchronize(stream);
    }

    // apply StandardWells
    if (this->num_std_wells > 0) {
        int smem_size = 2 * sizeof(Scalar) * this->dim_wells;
        apply_well_contributions <<< this->num_std_wells, 32, smem_size, stream>>>(d_Cnnzs,
                                                                                   d_Dnnzs,
                                                                                   d_Bnnzs,
                                                                                   d_Ccols,
                                                                                   d_Bcols,
                                                                                   d_x,
                                                                                   d_y,
                                                                                   this->dim,
                                                                                   this->dim_wells,
                                                                                   this->d_val_pointers);
    }
}

template<class Scalar>
void WellContributionsCuda<Scalar>::APIaddMatrix(MatrixType type, int* colIndices,
                                                 Scalar* values, unsigned int val_size)
{
    switch (type) {
    case MatrixType::C:
        cudaMemcpy(d_Cnnzs + this->num_blocks_so_far * this->dim * this->dim_wells,
                   values, sizeof(Scalar) * val_size * this->dim * this->dim_wells,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_Ccols + this->num_blocks_so_far, colIndices,
                   sizeof(int) * val_size, cudaMemcpyHostToDevice);
        break;
    case MatrixType::D:
        cudaMemcpy(d_Dnnzs + this->num_std_wells_so_far * this->dim_wells * this->dim_wells,
                   values, sizeof(Scalar) * this->dim_wells * this->dim_wells,
                   cudaMemcpyHostToDevice);
        break;
    case MatrixType::B:
        cudaMemcpy(d_Bnnzs + this->num_blocks_so_far * this->dim * this->dim_wells,
                   values, sizeof(Scalar) * val_size * this->dim * this->dim_wells,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(d_Bcols + this->num_blocks_so_far, colIndices,
                   sizeof(int) * val_size, cudaMemcpyHostToDevice);
        this->val_pointers[this->num_std_wells_so_far] = this->num_blocks_so_far;
        if (this->num_std_wells_so_far == this->num_std_wells - 1) {
            this->val_pointers[this->num_std_wells] = this->num_blocks;
            cudaMemcpy(d_val_pointers, this->val_pointers.data(),
                       sizeof(unsigned int) * (this->num_std_wells + 1),
                       cudaMemcpyHostToDevice);
        }
        break;
    default:
        OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributions::addMatrix()");
    }
    cudaCheckLastError("WellContributions::addMatrix() failed");
}

template<class Scalar>
void WellContributionsCuda<Scalar>::setCudaStream(cudaStream_t stream_)
{
    this->stream = stream_;
    for (auto& well : this->multisegments) {
        well->setCudaStream(stream_);
    }
}

template class WellContributionsCuda<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellContributionsCuda<float>;
#endif

} //namespace Opm

