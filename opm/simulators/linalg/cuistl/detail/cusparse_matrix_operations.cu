/*
  Copyright 2022-2023 SINTEF AS

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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <stdexcept>
#include <config.h>

namespace Opm::cuistl::detail
{
namespace
{

    // TODO: figure out if this can be generalized effectively, this seems excessively verbose
    // explicit formulas based on Dune cpu code
    template <class T, int blocksize>
    __device__ __forceinline__ void invBlockOutOfPlace(const T* __restrict__ srcBlock, T* __restrict__ dstBlock)
    {
        if (blocksize == 1) {
            dstBlock[0] = 1.0 / (srcBlock[0]);
        } else if (blocksize == 2) {
            T detInv = 1.0 / (srcBlock[0] * srcBlock[3] - srcBlock[1] * srcBlock[2]);
            dstBlock[0] = srcBlock[3] * detInv;
            dstBlock[1] = -srcBlock[1] * detInv;
            dstBlock[2] = -srcBlock[2] * detInv;
            dstBlock[3] = srcBlock[0] * detInv;
        } else if (blocksize == 3) {
            // based on Dune implementation
            T t4 = srcBlock[0] * srcBlock[4];
            T t6 = srcBlock[0] * srcBlock[5];
            T t8 = srcBlock[1] * srcBlock[3];
            T t10 = srcBlock[2] * srcBlock[3];
            T t12 = srcBlock[1] * srcBlock[6];
            T t14 = srcBlock[2] * srcBlock[6];

            T t17 = 1.0
                / (t4 * srcBlock[8] - t6 * srcBlock[7] - t8 * srcBlock[8] + t10 * srcBlock[7] + t12 * srcBlock[5]
                   - t14 * srcBlock[4]); // t17 is 1/determinant

            dstBlock[0] = (srcBlock[4] * srcBlock[8] - srcBlock[5] * srcBlock[7]) * t17;
            dstBlock[1] = -(srcBlock[1] * srcBlock[8] - srcBlock[2] * srcBlock[7]) * t17;
            dstBlock[2] = (srcBlock[1] * srcBlock[5] - srcBlock[2] * srcBlock[4]) * t17;
            dstBlock[3] = -(srcBlock[3] * srcBlock[8] - srcBlock[5] * srcBlock[6]) * t17;
            dstBlock[4] = (srcBlock[0] * srcBlock[8] - t14) * t17;
            dstBlock[5] = -(t6 - t10) * t17;
            dstBlock[6] = (srcBlock[3] * srcBlock[7] - srcBlock[4] * srcBlock[6]) * t17;
            dstBlock[7] = -(srcBlock[0] * srcBlock[7] - t12) * t17;
            dstBlock[8] = (t4 - t8) * t17;
        }
    }

    // explicit formulas based on Dune cpu code
    template <class T, int blocksize>
    __device__ __forceinline__ void invBlockInPlace(T* __restrict__ block)
    {
        if (blocksize == 1) {
            block[0] = 1.0 / (block[0]);
        } else if (blocksize == 2) {
            T detInv = 1.0 / (block[0] * block[3] - block[1] * block[2]);

            T temp = block[0];
            block[0] = block[3] * detInv;
            block[1] = -block[1] * detInv;
            block[2] = -block[2] * detInv;
            block[3] = temp * detInv;
        } else if (blocksize == 3) {
            T t4 = block[0] * block[4];
            T t6 = block[0] * block[5];
            T t8 = block[1] * block[3];
            T t10 = block[2] * block[3];
            T t12 = block[1] * block[6];
            T t14 = block[2] * block[6];

            T det = (t4 * block[8] - t6 * block[7] - t8 * block[8] + t10 * block[7] + t12 * block[5] - t14 * block[4]);
            T t17 = T(1.0) / det;

            T matrix01 = block[1];
            T matrix00 = block[0];
            T matrix10 = block[3];
            T matrix11 = block[4];

            block[0] = (block[4] * block[8] - block[5] * block[7]) * t17;
            block[1] = -(block[1] * block[8] - block[2] * block[7]) * t17;
            block[2] = (matrix01 * block[5] - block[2] * block[4]) * t17;
            block[3] = -(block[3] * block[8] - block[5] * block[6]) * t17;
            block[4] = (matrix00 * block[8] - t14) * t17;
            block[5] = -(t6 - t10) * t17;
            block[6] = (matrix10 * block[7] - matrix11 * block[6]) * t17;
            block[7] = -(matrix00 * block[7] - t12) * t17;
            block[8] = (t4 - t8) * t17;
        }
    }

    enum class MVType { SET, PLUS, MINUS };
    // SET:   c  = A*b
    // PLS:   c += A*b
    // MINUS: c -= A*b
    template <class T, int blocksize, MVType OP>
    __device__ __forceinline__ void matrixVectorProductWithAction(const T* A, const T* b, T* c)
    {
        for (int i = 0; i < blocksize; ++i) {
            if (OP == MVType::SET) {
                c[i] = 0;
            }

            for (int j = 0; j < blocksize; ++j) {
                if (OP == MVType::SET || OP == MVType::PLUS) {
                    c[i] += A[i * blocksize + j] * b[j];
                } else if (OP == MVType::MINUS) {
                    c[i] -= A[i * blocksize + j] * b[j];
                }
            }
        }
    }

    template <class T, int blocksize>
    __device__ __forceinline__ void mv(const T* a, const T* b, T* c)
    {
        matrixVectorProductWithAction<T, blocksize, MVType::SET>(a, b, c);
    }

    template <class T, int blocksize>
    __device__ __forceinline__ void umv(const T* a, const T* b, T* c)
    {
        matrixVectorProductWithAction<T, blocksize, MVType::PLUS>(a, b, c);
    }

    template <class T, int blocksize>
    __device__ __forceinline__ void mmv(const T* a, const T* b, T* c)
    {
        matrixVectorProductWithAction<T, blocksize, MVType::MINUS>(a, b, c);
    }

    // dst -= A*B*C
    template <class T, int blocksize>
    __device__ __forceinline__ void mmx2Subtraction(T* A, T* B, T* C, T* dst)
    {

        T tmp[blocksize * blocksize] = {0};
        // tmp = A*B
        for (int i = 0; i < blocksize; ++i) {
            for (int k = 0; k < blocksize; ++k) {
                for (int j = 0; j < blocksize; ++j) {
                    tmp[i * blocksize + j] += A[i * blocksize + k] * B[k * blocksize + j];
                }
            }
        }

        // dst = tmp*C
        for (int i = 0; i < blocksize; ++i) {
            for (int k = 0; k < blocksize; ++k) {
                for (int j = 0; j < blocksize; ++j) {
                    dst[i * blocksize + j] -= tmp[i * blocksize + k] * C[k * blocksize + j];
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void
    cuInvertDiagonalAndFlatten(T* matNonZeroValues, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;

        if (row < numberOfRows) {
            size_t nnzIdx = rowIndices[row];
            size_t nnzIdxLim = rowIndices[row + 1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a
            // zero element
            while (colIndices[nnzIdx] != row && nnzIdx <= nnzIdxLim) {
                ++nnzIdx;
            }

            // diagBlock points to the start of where the diagonal block is stored
            T* diagBlock = &matNonZeroValues[blocksize * blocksize * nnzIdx];
            // vecBlock points to the start of the block element in the vector where the inverse of the diagonal block
            // element should be stored
            T* vecBlock = &vec[blocksize * blocksize * row];

            invBlockOutOfPlace<T, blocksize>(diagBlock, vecBlock);
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeLowerSolveLevelSet(T* mat,
                                                int* rowIndices,
                                                int* colIndices,
                                                int* indexConversion,
                                                int startIdx,
                                                int rowsInLevelSet,
                                                const T* dInv,
                                                const T* d,
                                                T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; colIndices[block] < naturalRowIdx; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeLowerSolveLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     const T* d,
                                                     T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            // TODO: removce the first condition in the for loop
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }


    template <class T, int blocksize>
    __global__ void cuComputeUpperSolveLevelSet(T* mat,
                                                int* rowIndices,
                                                int* colIndices,
                                                int* indexConversion,
                                                int startIdx,
                                                int rowsInLevelSet,
                                                const T* dInv,
                                                T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize] = {0};
            for (int block = nnzIdxLim - 1; colIndices[block] > naturalRowIdx; --block) {
                const int col = colIndices[block];
                umv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mmv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeUpperSolveLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize] = {0};
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                umv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mmv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeDiluDiagonal(T* mat,
                                          int* rowIndices,
                                          int* colIndices,
                                          int* reorderedToNatural,
                                          int* naturalToReordered,
                                          const int startIdx,
                                          int rowsInLevelSet,
                                          T* dInv)
    {
        const auto reorderedRowIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const int naturalRowIdx = reorderedToNatural[reorderedRowIdx];
            const size_t nnzIdx = rowIndices[reorderedRowIdx];

            int diagIdx = nnzIdx;
            while (colIndices[diagIdx] != naturalRowIdx) {
                ++diagIdx;
            }

            T dInvTmp[blocksize * blocksize];
            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInvTmp[i * blocksize + j] = mat[diagIdx * blocksize * blocksize + i * blocksize + j];
                }
            }

            for (int block = nnzIdx; colIndices[block] < naturalRowIdx; ++block) {
                const int col = naturalToReordered[colIndices[block]];
                // find element with indices swapped
                // Binary search over block in the right row, [rowIndices[col], rowindices[col+1]-1] defines the range
                // we binary search over
                int left = rowIndices[col];
                int right = rowIndices[col + 1] - 1;
                int mid = left + (right - left) / 2; // overflow-safe average;

                while (left <= right) {
                    const int col = colIndices[mid];

                    if (col == naturalRowIdx) {
                        break;
                    } else if (col < naturalRowIdx) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                    mid = left + (right - left) / 2; // overflow-safe average
                }

                const int symOpposite = mid;

                mmx2Subtraction<T, blocksize>(&mat[block * blocksize * blocksize],
                                              &dInv[col * blocksize * blocksize],
                                              &mat[symOpposite * blocksize * blocksize],
                                              dInvTmp);
            }

            invBlockInPlace<T, blocksize>(dInvTmp);

            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInv[reorderedRowIdx * blocksize * blocksize + i * blocksize + j] = dInvTmp[i * blocksize + j];
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeDiluDiagonalSplit(T* reorderedLowerMat,
                                               int* lowerRowIndices,
                                               int* lowerColIndices,
                                               T* reorderedUpperMat,
                                               int* upperRowIndices,
                                               int* upperColIndices,
                                               T* diagonal,
                                               int* reorderedToNatural,
                                               int* naturalToReordered,
                                               const int startIdx,
                                               int rowsInLevelSet,
                                               T* dInv)
    {
        const auto reorderedRowIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const int naturalRowIdx = reorderedToNatural[reorderedRowIdx];
            const size_t lowerRowStart = lowerRowIndices[reorderedRowIdx];
            const size_t lowerRowEnd = lowerRowIndices[reorderedRowIdx + 1];

            T dInvTmp[blocksize * blocksize];
            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInvTmp[i * blocksize + j] = diagonal[reorderedRowIdx * blocksize * blocksize + i * blocksize + j];
                }
            }

            for (int block = lowerRowStart; block < lowerRowEnd; ++block) {
                const int col = naturalToReordered[lowerColIndices[block]];

                int symOppositeIdx = upperRowIndices[col];
                for (; symOppositeIdx < upperRowIndices[col + 1]; ++symOppositeIdx) {
                    if (naturalRowIdx == upperColIndices[symOppositeIdx]) {
                        break;
                    }
                }

                const int symOppositeBlock = symOppositeIdx;

                mmx2Subtraction<T, blocksize>(&reorderedLowerMat[block * blocksize * blocksize],
                                              &dInv[col * blocksize * blocksize],
                                              &reorderedUpperMat[symOppositeBlock * blocksize * blocksize],
                                              dInvTmp);
            }

            invBlockInPlace<T, blocksize>(dInvTmp);

            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInv[reorderedRowIdx * blocksize * blocksize + i * blocksize + j] = dInvTmp[i * blocksize + j];
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuMoveDataToReordered(
        T* srcMatrix, int* srcRowIndices, T* dstMatrix, int* dstRowIndices, int* indexConversion, size_t numberOfRows)
    {
        const auto srcRow = blockDim.x * blockIdx.x + threadIdx.x;
        if (srcRow < numberOfRows) {

            const auto dstRow = indexConversion[srcRow];

            for (int srcBlock = srcRowIndices[srcRow], dstBlock = dstRowIndices[dstRow];
                 srcBlock < srcRowIndices[srcRow + 1];
                 ++srcBlock, ++dstBlock) {
                for (int i = 0; i < blocksize; ++i) {
                    for (int j = 0; j < blocksize; ++j) {
                        dstMatrix[dstBlock * blocksize * blocksize + i * blocksize + j]
                            = srcMatrix[srcBlock * blocksize * blocksize + i * blocksize + j];
                    }
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuMoveDataToReorderedSplit(T* srcMatrix,
                                               int* srcRowIndices,
                                               int* srcColumnIndices,
                                               T* dstLowerMatrix,
                                               int* dstLowerRowIndices,
                                               T* dstUpperMatrix,
                                               int* dstUpperRowIndices,
                                               T* dstDiag,
                                               int* naturalToReordered,
                                               size_t numberOfRows)
    {
        const auto srcRow = blockDim.x * blockIdx.x + threadIdx.x;
        if (srcRow < numberOfRows) {

            const auto dstRow = naturalToReordered[srcRow];
            const auto rowStart = srcRowIndices[srcRow];
            const auto rowEnd = srcRowIndices[srcRow + 1];

            auto lowerBlock = dstLowerRowIndices[dstRow];
            auto upperBlock = dstUpperRowIndices[dstRow];

            for (int srcBlock = rowStart; srcBlock < rowEnd; srcBlock++) {
                int dstBlock;
                T* dstBuffer;

                if (srcColumnIndices[srcBlock] < srcRow) { // we are writing a value to the lower triangular matrix
                    dstBlock = lowerBlock;
                    ++lowerBlock;
                    dstBuffer = dstLowerMatrix;
                } else if (srcColumnIndices[srcBlock]
                           > srcRow) { // we are writing a value to the upper triangular matrix
                    dstBlock = upperBlock;
                    ++upperBlock;
                    dstBuffer = dstUpperMatrix;
                } else { // we are writing a value to the diagonal
                    dstBlock = dstRow;
                    dstBuffer = dstDiag;
                }
                for (int i = 0; i < blocksize; ++i) {
                    for (int j = 0; j < blocksize; ++j) {
                        dstBuffer[dstBlock * blocksize * blocksize + i * blocksize + j]
                            = srcMatrix[srcBlock * blocksize * blocksize + i * blocksize + j];
                    }
                }
            }
        }
    }

    // Kernel here is the function object of the cuda kernel
    template <class Kernel>
    inline int getCudaRecomendedThreadBlockSize(Kernel k, int suggestedThrBlockSize=-1)
    {
        if (suggestedThrBlockSize != -1){
            return suggestedThrBlockSize;
        }
// Use cuda API to maximize occupancy, otherwise we just pick a thread block size if it is not tuned
#if USE_HIP
    return 512;
#else
        int blockSize = 1024;
        int tmpGridSize;
        std::ignore = cudaOccupancyMaxPotentialBlockSize(&tmpGridSize, &blockSize, k, 0, 0);
        return blockSize;
#endif
    }

    inline int getNumberOfBlocks(int wantedThreads, int threadBlockSize)
    {
        return (wantedThreads + threadBlockSize - 1) / threadBlockSize;
    }

} // namespace

template <class T, int blocksize>
void
invertDiagonalAndFlatten(T* mat, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec)
{
    if (blocksize <= 3) {
        int threadBlockSize = getCudaRecomendedThreadBlockSize(cuInvertDiagonalAndFlatten<T, blocksize>);
        int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);
        cuInvertDiagonalAndFlatten<T, blocksize>
            <<<nThreadBlocks, threadBlockSize>>>(mat, rowIndices, colIndices, numberOfRows, vec);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

// perform the lower solve for all rows in the same level set
template <class T, int blocksize>
void
computeLowerSolveLevelSet(T* reorderedMat,
                          int* rowIndices,
                          int* colIndices,
                          int* indexConversion,
                          int startIdx,
                          int rowsInLevelSet,
                          const T* dInv,
                          const T* d,
                          T* v,
                          int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeLowerSolveLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuComputeLowerSolveLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}


template <class T, int blocksize>
void
computeLowerSolveLevelSetSplit(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               const T* d,
                               T* v,
                               int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeLowerSolveLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuComputeLowerSolveLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}
// perform the upper solve for all rows in the same level set
template <class T, int blocksize>
void
computeUpperSolveLevelSet(T* reorderedMat,
                          int* rowIndices,
                          int* colIndices,
                          int* indexConversion,
                          int startIdx,
                          int rowsInLevelSet,
                          const T* dInv,
                          T* v,
                          int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeUpperSolveLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuComputeUpperSolveLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <class T, int blocksize>
void
computeUpperSolveLevelSetSplit(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               T* v,
                               int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeUpperSolveLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuComputeUpperSolveLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <class T, int blocksize>
void
computeDiluDiagonal(T* reorderedMat,
                    int* rowIndices,
                    int* colIndices,
                    int* reorderedToNatural,
                    int* naturalToReordered,
                    const int startIdx,
                    int rowsInLevelSet,
                    T* dInv,
                    int thrBlockSize)
{
    if (blocksize <= 3) {
        int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeDiluDiagonal<T, blocksize>, thrBlockSize);
        int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        cuComputeDiluDiagonal<T, blocksize>
            <<<nThreadBlocks, threadBlockSize>>>(reorderedMat,
                                                                        rowIndices,
                                                                        colIndices,
                                                                        reorderedToNatural,
                                                                        naturalToReordered,
                                                                        startIdx,
                                                                        rowsInLevelSet,
                                                                        dInv);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

template <class T, int blocksize>
void
computeDiluDiagonalSplit(T* reorderedLowerMat,
                         int* lowerRowIndices,
                         int* lowerColIndices,
                         T* reorderedUpperMat,
                         int* upperRowIndices,
                         int* upperColIndices,
                         T* diagonal,
                         int* reorderedToNatural,
                         int* naturalToReordered,
                         const int startIdx,
                         int rowsInLevelSet,
                         T* dInv,
                         int thrBlockSize)
{
    if (blocksize <= 3) {
        int threadBlockSize = getCudaRecomendedThreadBlockSize(cuComputeDiluDiagonalSplit<T, blocksize>, thrBlockSize);
        int nThreadBlocks = getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        cuComputeDiluDiagonalSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(reorderedLowerMat,
                                                                                     lowerRowIndices,
                                                                                     lowerColIndices,
                                                                                     reorderedUpperMat,
                                                                                     upperRowIndices,
                                                                                     upperColIndices,
                                                                                     diagonal,
                                                                                     reorderedToNatural,
                                                                                     naturalToReordered,
                                                                                     startIdx,
                                                                                     rowsInLevelSet,
                                                                                     dInv);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

template <class T, int blocksize>
void
copyMatDataToReordered(
    T* srcMatrix, int* srcRowIndices, T* dstMatrix, int* dstRowIndices, int* naturalToReordered, size_t numberOfRows,
    int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuMoveDataToReordered<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);
    cuMoveDataToReordered<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        srcMatrix, srcRowIndices, dstMatrix, dstRowIndices, naturalToReordered, numberOfRows);
}

template <class T, int blocksize>
void
copyMatDataToReorderedSplit(T* srcMatrix,
                            int* srcRowIndices,
                            int* srcColumnIndices,
                            T* dstLowerMatrix,
                            int* dstLowerRowIndices,
                            T* dstUpperMatrix,
                            int* dstUpperRowIndices,
                            T* dstDiag,
                            int* naturalToReordered,
                            size_t numberOfRows,
                            int thrBlockSize)
{
    int threadBlockSize = getCudaRecomendedThreadBlockSize(cuMoveDataToReorderedSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);
    cuMoveDataToReorderedSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(srcMatrix,
                                                                                 srcRowIndices,
                                                                                 srcColumnIndices,
                                                                                 dstLowerMatrix,
                                                                                 dstLowerRowIndices,
                                                                                 dstUpperMatrix,
                                                                                 dstUpperRowIndices,
                                                                                 dstDiag,
                                                                                 naturalToReordered,
                                                                                 numberOfRows);
}

#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void invertDiagonalAndFlatten<T, blocksize>(T*, int*, int*, size_t, T*);                                  \
    template void copyMatDataToReordered<T, blocksize>(T*, int*, T*, int*, int*, size_t, int);                              \
    template void copyMatDataToReorderedSplit<T, blocksize>(T*, int*, int*, T*, int*, T*, int*, T*, int*, size_t, int);     \
    template void computeDiluDiagonal<T, blocksize>(T*, int*, int*, int*, int*, const int, int, T*, int);                   \
    template void computeDiluDiagonalSplit<T, blocksize>(                                                              \
        T*, int*, int*, T*, int*, int*, T*, int*, int*, const int, int, T*, int);                                           \
    template void computeUpperSolveLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);               \
    template void computeLowerSolveLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, const T*, T*, int);     \
    template void computeUpperSolveLevelSetSplit<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);          \
    template void computeLowerSolveLevelSetSplit<T, blocksize>(T*, int*, int*, int*, int, int, const T*, const T*, T*, int);

INSTANTIATE_KERNEL_WRAPPERS(float, 1);
INSTANTIATE_KERNEL_WRAPPERS(float, 2);
INSTANTIATE_KERNEL_WRAPPERS(float, 3);
INSTANTIATE_KERNEL_WRAPPERS(float, 4);
INSTANTIATE_KERNEL_WRAPPERS(float, 5);
INSTANTIATE_KERNEL_WRAPPERS(float, 6);
INSTANTIATE_KERNEL_WRAPPERS(double, 1);
INSTANTIATE_KERNEL_WRAPPERS(double, 2);
INSTANTIATE_KERNEL_WRAPPERS(double, 3);
INSTANTIATE_KERNEL_WRAPPERS(double, 4);
INSTANTIATE_KERNEL_WRAPPERS(double, 5);
INSTANTIATE_KERNEL_WRAPPERS(double, 6);
} // namespace Opm::cuistl::detail
