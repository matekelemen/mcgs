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
    // Allocate arrays for the new partition
    std::vector<TIndex> newPartitionExtents(pPartition->size() + 1);
    newPartitionExtents[0] = static_cast<TIndex>(0);

    // Allocate temporary matrix
    const auto partitionCount = pPartition->size();
    std::vector<TIndex> newRowExtents(rowCount + 1);
    std::vector<TIndex> newColumnIndices(nonzeroCount);
    std::vector<TValue> newNonzeros(nonzeroCount);
    std::vector<TValue> rhs(rowCount);
    newRowExtents.front() = 0;

    // Compute new extents
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

            newPartitionExtents[iPartition + 1] = newPartitionExtents[iPartition] + partitionSize;
        } // for iPartition in range(partitionCount)
    }

    std::vector<TIndex> columnMap(columnCount);

    #ifdef MCGS_OPENMP
    #pragma omp parallel
    #endif
    {
        for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
            auto itPartitionBegin = pPartition->begin(iPartition);
            const auto partitionSize = pPartition->size(iPartition);

            #ifdef MCGS_OPENMP
            #pragma omp for
            #endif
            for (std::remove_const_t<decltype(partitionSize)> iLocal=0; iLocal<partitionSize; ++iLocal) {
                const TIndex iOldRow = itPartitionBegin[iLocal];
                const TIndex iNewRow = newPartitionExtents[iPartition] + iLocal;
                columnMap[iOldRow] = iNewRow;

                std::copy(pColumnIndices + pRowExtents[iOldRow],
                          pColumnIndices + pRowExtents[iOldRow + 1],
                          newColumnIndices.data() + newRowExtents[iNewRow]);

                std::copy(pNonzeros + pRowExtents[iOldRow],
                          pNonzeros + pRowExtents[iOldRow + 1],
                          newNonzeros.data() + newRowExtents[iNewRow]);

                rhs[iNewRow] = pRHS[iOldRow];
            } // for iLocal in range(parititionSize)
        } // for iPartition in range(partitionCount)
    } // omp parallel

    std::copy(newRowExtents.begin(), newRowExtents.end(), pRowExtents);
    std::copy(newNonzeros.begin(), newNonzeros.end(), pNonzeros);
    std::copy(rhs.begin(), rhs.end(), pRHS);
    std::transform(newColumnIndices.begin(),
                   newColumnIndices.end(),
                   pColumnIndices,
                   [&columnMap](const TIndex iOldColumn) -> TIndex {
                        return columnMap[iOldColumn];
                   });

    // Allocate arrays for the new partition
    // and assign the new row indices
    std::vector<TIndex> newRowIndices(rowCount + 1);
    std::iota(newRowIndices.begin(), newRowIndices.end(), TIndex(0));

    return new Partition<TIndex>(std::move(newPartitionExtents), std::move(newRowIndices));
}


template <class TIndex, class TValue>
int revertReorder(TValue* pRHS,
                  const TIndex columnCount,
                  const Partition<TIndex>* pPartition)
{
    std::vector<TValue> swap(columnCount);

    TIndex iNewRow = 0;
    for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
        for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            swap[*itPartition] = pRHS[iNewRow++];
        }
    }

    if (iNewRow != columnCount) {
        return MCGS_FAILURE;
    }

    std::copy(swap.begin(), swap.end(), pRHS);
    return MCGS_SUCCESS;
}


template <class TIndex, class TValue>
int revertReorder(const TIndex rowCount, const TIndex columnCount, const TIndex nonzeroCount,
                  TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
                  TValue* pRHS,
                  const Partition<TIndex>* pPartition)
{
    // Construct inverse partition
    std::vector<TIndex> partitionExtents {static_cast<TIndex>(0), rowCount}, rowIndices(rowCount);
    const TIndex* itPartitionBegin = pPartition->begin(0);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (TIndex iRow=0; iRow<rowCount; ++iRow) {
        rowIndices[itPartitionBegin[iRow]] = iRow;
    } // for iRow in range(rowCount)

    Partition<TIndex> inversePartition(std::move(partitionExtents), std::move(rowIndices));

    // Reorder
    Partition<TIndex>* pDummy = reorder(rowCount, columnCount, nonzeroCount,
                                        pRowExtents, pColumnIndices, pNonzeros,
                                        pRHS,
                                        &inversePartition);

    if (pDummy == nullptr) {
        return MCGS_FAILURE;
    } else {
        delete pDummy;
    }

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_REORDER(TIndex, TValue)                                                    \
    template Partition<TIndex>* reorder<TIndex,TValue>(const TIndex, const TIndex, const TIndex,    \
                                        TIndex*, TIndex*, TValue*,                                  \
                                        TValue*,                                                    \
                                        const Partition<TIndex>*);                                  \
    template int revertReorder<TIndex,TValue>(TValue*,                                              \
                                              const TIndex,                                         \
                                              const Partition<TIndex>*);                            \
    template int revertReorder<TIndex,TValue>(const TIndex, const TIndex, const TIndex,             \
                                              TIndex*, TIndex*, TValue*,                            \
                                              TValue*,                                              \
                                              const Partition<TIndex>*)

MCGS_INSTANTIATE_REORDER(int, double);
MCGS_INSTANTIATE_REORDER(long, double);
MCGS_INSTANTIATE_REORDER(unsigned, double);
MCGS_INSTANTIATE_REORDER(std::size_t, double);

#undef MCGS_INSTANTIATE_REORDER


} // namespace mcgs
