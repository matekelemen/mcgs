// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::partition, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <cstddef> // std::size_t


namespace mcgs {


template <class TIndex, class TColor>
Partition<TIndex,TColor>::Partition(const TColor* pColor, const TIndex columnCount) noexcept
    : _partitionExtents(1, 0),
      _partitions()
{
    // @todo
}


template <class TIndex, class TColor>
typename Partition<TIndex,TColor>::size_type
Partition<TIndex,TColor>::size() const noexcept
{
    return _partitionExtents.size() - 1;
}


template <class TIndex, class TColor>
typename Partition<TIndex,TColor>::size_type
Partition<TIndex,TColor>::size(const size_type iPartition) const noexcept
{
    return std::distance(this->begin(iPartition), this->end(iPartition));
}


template <class TIndex, class TColor>
typename Partition<TIndex,TColor>::const_iterator
Partition<TIndex,TColor>::begin(const size_type iPartition) const noexcept
{
    return &_partitions[_partitionExtents[iPartition]];
}


template <class TIndex, class TColor>
typename Partition<TIndex,TColor>::const_iterator
Partition<TIndex,TColor>::end(const size_type iPartition) const noexcept
{
    return &_partitions[_partitionExtents[iPartition + 1]];
}


template <class TIndex, class TColor>
[[nodiscard]] Partition<TIndex,TColor>* makePartition(const TColor* pColors,
                                                      const TIndex columnCount)
{
    return new Partition<TIndex,TColor>(pColors, columnCount);
}


template <class TIndex, class TColor>
void destroyPartition(Partition<TIndex,TColor>* pPartition)
{
    delete pPartition;
}


#define MCGS_INSTANTIATE_PARTITION(TIndex, TColor)                                                  \
    template Partition<TIndex,TColor>* makePartition<TIndex,TColor>(const TColor*, const TIndex);   \
    template void destroyPartition<TIndex,TColor>(Partition<TIndex,TColor>*)

MCGS_INSTANTIATE_PARTITION(int, unsigned);

MCGS_INSTANTIATE_PARTITION(long, unsigned);

MCGS_INSTANTIATE_PARTITION(unsigned, unsigned);

MCGS_INSTANTIATE_PARTITION(std::size_t, unsigned);

MCGS_INSTANTIATE_PARTITION(std::size_t, std::size_t);

#undef MCGS_INSTANTIATE_PARTITION


} // namespace mcgs
