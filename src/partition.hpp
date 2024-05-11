#pragma once

// --- STL Includes ---
#include <vector> // std::vector
#include <cstddef> // std::size_t, std::ptrdiff_t


namespace mcgs {


template <class TIndex, class TColor>
class Partition
{
public:
    using value_type = TIndex;
    using reference = TIndex&;
    using const_reference = const TIndex&;
    using iterator = const TIndex*;
    using const_iterator = const TIndex*;
    using difference_type = std::ptrdiff_t;
    using size_type = std::size_t;

    Partition(const TColor* pColors, const TIndex columnCount) noexcept;

    size_type size() const noexcept;

    size_type size(const size_type iPartition) const noexcept;

    const_iterator begin(const size_type iPartition) const noexcept;

    const_iterator end(const size_type iPartition) const noexcept;

private:
    std::vector<TIndex> _partitionExtents;

    std::vector<TIndex> _partitions;

    Partition() = delete;
    Partition(Partition&&) = delete;
    Partition(const Partition&) = delete;
    Partition& operator=(Partition&&) = delete;
    Partition& operator=(const Partition&) = delete;
};


} // namespace mcgs
