// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::reorder, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <vector> // std::vector
#include <algorithm> // std::copy, std::swap
#include <numeric> // std::iota


namespace mcgs {


template <class TIndex, class TValue, class TColor>
int reorder(const TIndex rowCount, const TIndex columnCount, const TIndex nonzeroCount,
            TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
            TValue* pRHS,
            Partition<TIndex,TColor>* pPartition)
{
    const auto partitionCount = pPartition->size();
    std::vector<TIndex> newRowExtents(rowCount + 1);
    std::vector<TIndex> newColumnIndices(nonzeroCount);
    std::vector<TValue> newNonzeros(nonzeroCount);
    newRowExtents.front() = 0;

    {
        TIndex iNewRow = 0;

        for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
            const auto itPartitionBegin = pPartition->begin(iPartition);
            const auto partitionSize = pPartition->size(iPartition);

            for (std::remove_const_t<decltype(partitionSize)> iLocal=0; iLocal<partitionSize; ++iLocal) {
                const TIndex iOldRow = itPartitionBegin[iLocal];
                const auto rowSize = pRowExtents[iOldRow + 1] - pRowExtents[iOldRow];
                newRowExtents[iNewRow + 1] = newRowExtents[iNewRow] + rowSize;
                ++iNewRow;
            } // for iLocal in range(parititionSize)
        } // for iPartition in range(partitionCount)
    }

    TIndex iNewRowBegin = 0;
    for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
        auto itPartitionBegin = pPartition->begin(iPartition);
        const auto partitionSize = pPartition->size(iPartition);

        //#pragma omp parallel for
        for (std::remove_const_t<decltype(partitionSize)> iLocal=0; iLocal<partitionSize; ++iLocal) {
            const TIndex iOldRow = itPartitionBegin[iLocal];
            const TIndex iNewRow = iNewRowBegin + iLocal;

            std::copy(pColumnIndices + pRowExtents[iOldRow],
                      pColumnIndices + pRowExtents[iOldRow + 1],
                      newColumnIndices.data() + newRowExtents[iNewRow]);

            std::copy(pNonzeros + pRowExtents[iOldRow],
                      pNonzeros + pRowExtents[iOldRow + 1],
                      newNonzeros.data() + newRowExtents[iNewRow]);

            std::swap(pRHS[iOldRow], pRHS[iNewRow]);
        } // for iLocal in range(parititionSize)

        std::iota(itPartitionBegin, itPartitionBegin + partitionSize, iNewRowBegin);
        iNewRowBegin += partitionSize;
    } // for iPartition in range(partitionCount)

    std::copy(newRowExtents.begin(), newRowExtents.end(), pRowExtents);
    std::copy(newColumnIndices.begin(), newColumnIndices.end(), pColumnIndices);
    std::copy(newNonzeros.begin(), newNonzeros.end(), pNonzeros);

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_REORDER(TIndex, TValue, TColor)                                    \
    template int reorder<TIndex,TValue,TColor>(const TIndex, const TIndex, const TIndex,    \
                                               TIndex*, TIndex*, TValue*,                   \
                                               TValue*,                                     \
                                               Partition<TIndex,TColor>*);

MCGS_INSTANTIATE_REORDER(int, double, unsigned);

MCGS_INSTANTIATE_REORDER(long, double, unsigned);

MCGS_INSTANTIATE_REORDER(unsigned, double, unsigned);

MCGS_INSTANTIATE_REORDER(std::size_t, double, unsigned);

MCGS_INSTANTIATE_REORDER(std::size_t, double, std::size_t);

#undef MCGS_INSTANTIATE_REORDER


} // namespace mcgs
