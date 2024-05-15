// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::partition, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <cstddef> // std::size_t
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <algorithm> // std::copy


namespace mcgs {


template <class TIndex>
template <class TColor>
Partition<TIndex>::Partition(const TColor* pColors, const TIndex columnCount)
    : _partitionExtents(1, 0),
      _rowIndices()
{
    std::unordered_map<
        TColor,
        std::vector<TIndex>
    > partitions;

    for (TIndex iColumn=0; iColumn<columnCount; ++iColumn) {
        const TColor color = pColors[iColumn];
        partitions.emplace(color, std::vector<TIndex> {})   // <== make sure an entry is mapped to color
            .first                                          // <== iterator pointing to the entry
            ->second                                        // <== reference to the mapped vector
            .push_back(iColumn);                            // <== insert the column index into the mapped vector
    } // for iColumn in range(columnCount)

    _rowIndices.reserve(columnCount);
    for ([[maybe_unused]] const auto& [iColor, rColumns] : partitions) {
        const std::size_t iPartitionBegin = _rowIndices.size();
        _rowIndices.resize(iPartitionBegin + rColumns.size());

        _partitionExtents.push_back(_rowIndices.size());
        std::copy(rColumns.begin(), rColumns.end(), _rowIndices.begin() + iPartitionBegin);
    }
}


template <class TIndex>
Partition<TIndex>::Partition(std::vector<TIndex>&& rPartitionExtents,
                             std::vector<TIndex>&& rRowIndices) noexcept
    : _partitionExtents(std::move(rPartitionExtents)),
      _rowIndices(std::move(rRowIndices))
{}


template <class TIndex, class TColor>
[[nodiscard]] Partition<TIndex>* makePartition(const TColor* pColors,
                                               const TIndex columnCount)
{
    return new Partition<TIndex>(pColors, columnCount);
}


template <class TIndex>
void destroyPartition(Partition<TIndex>* pPartition)
{
    delete pPartition;
}


#define MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, TColor)                                  \
    template Partition<TIndex>* makePartition<TIndex,TColor>(const TColor*, const TIndex);

#define MCGS_INSTANTIATE_PARTITION(TIndex)                      \
    MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, unsigned)        \
    MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, std::size_t)     \
    template void destroyPartition<TIndex>(Partition<TIndex>*); \
    template class Partition<TIndex>

MCGS_INSTANTIATE_PARTITION(int);
MCGS_INSTANTIATE_PARTITION(long);
MCGS_INSTANTIATE_PARTITION(unsigned);
MCGS_INSTANTIATE_PARTITION(std::size_t);

#undef MCGS_INSTANTIATE_PARTITION_FACTORY
#undef MCGS_INSTANTIATE_PARTITION


} // namespace mcgs
