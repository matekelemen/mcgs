// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::reorder, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <vector> // std::vector
#include <algorithm> // std::copy, std::swap, std::transform
#include <numeric> // std::iota


namespace mcgs {


template <class TIndex, class TValue>
Partition<TIndex>* reorder(const TIndex rowCount, const TIndex columnCount, const TIndex nonzeroCount,
                           TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
                           TValue* pRHS,
                           const Partition<TIndex>* pPartition)
{
    // Allocate array for the new partition extents
    std::vector<TIndex> newPartitionExtents(pPartition->size() + 1);

    // Allocate temporary matrix
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
    std::vector<TIndex> columnMap(columnCount);

    for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
        auto itPartitionBegin = pPartition->begin(iPartition);
        const auto partitionSize = pPartition->size(iPartition);
        newPartitionExtents[iPartition] = iNewRowBegin;

        #pragma omp parallel for
        for (std::remove_const_t<decltype(partitionSize)> iLocal=0; iLocal<partitionSize; ++iLocal) {
            const TIndex iOldRow = itPartitionBegin[iLocal];
            const TIndex iNewRow = iNewRowBegin + iLocal;
            columnMap[iOldRow] = iNewRow;

            std::copy(pColumnIndices + pRowExtents[iOldRow],
                      pColumnIndices + pRowExtents[iOldRow + 1],
                      newColumnIndices.data() + newRowExtents[iNewRow]);

            std::copy(pNonzeros + pRowExtents[iOldRow],
                      pNonzeros + pRowExtents[iOldRow + 1],
                      newNonzeros.data() + newRowExtents[iNewRow]);

            std::swap(pRHS[iOldRow], pRHS[iNewRow]);
        } // for iLocal in range(parititionSize)

        iNewRowBegin += partitionSize;
    } // for iPartition in range(partitionCount)

    newPartitionExtents.back() = iNewRowBegin;
    std::copy(newRowExtents.begin(), newRowExtents.end(), pRowExtents);
    std::copy(newNonzeros.begin(), newNonzeros.end(), pNonzeros);
    std::transform(newColumnIndices.begin(),
                   newColumnIndices.end(),
                   pColumnIndices,
                   [&columnMap](const TIndex iOldColumn) -> TIndex {
                        return columnMap[iOldColumn];
                   });

    // Allocate arrays for the new partition
    // and assign the new row indices
    std::vector<TIndex> newRowIndices(rowCount);
    std::iota(newRowIndices.begin(), newRowIndices.end(), TIndex(0));

    return new Partition<TIndex>(std::move(newPartitionExtents), std::move(newRowIndices));
}


template <class TIndex, class TValue>
int revertReorder(TValue* pRHS, const Partition<TIndex>* pPartition)
{
    typename Partition<TIndex>::size_type iNewRow = 0;
    for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
        for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            std::swap(pRHS[iNewRow++], pRHS[*itPartition]);
        }
    }
    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_REORDER(TIndex, TValue)                                                    \
    template Partition<TIndex>* reorder<TIndex,TValue>(const TIndex, const TIndex, const TIndex,    \
                                        TIndex*, TIndex*, TValue*,                                  \
                                        TValue*,                                                    \
                                        const Partition<TIndex>*);                                  \
    template int revertReorder<TIndex,TValue>(TValue*,const Partition<TIndex>*)

MCGS_INSTANTIATE_REORDER(int, double);
MCGS_INSTANTIATE_REORDER(long, double);
MCGS_INSTANTIATE_REORDER(unsigned, double);
MCGS_INSTANTIATE_REORDER(std::size_t, double);

#undef MCGS_INSTANTIATE_REORDER


} // namespace mcgs
