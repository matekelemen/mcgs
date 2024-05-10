// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::GaussSeidel

// --- STL Includes ---
#include <vector> // std::vector
#include <set> // std::set, std::pair
#include <unordered_set> // std::unordered_set
#include <cstddef> // std::size_t
#include <algorithm> // std::min, std::max
#include <numeric> // std::iota
#include <random> // std::mt19937, std::uniform_int_distribution

#include <iostream>


namespace mcgs {


namespace detail {


template <class TIndex>
struct IndexPairTraits
{
    using value_type = std::pair<TIndex,TIndex>;

    struct Less
    {
        bool operator()(value_type left, value_type right) const noexcept
        {
            if (left.first < right.first) {
                return true;
            } else if (left.first == right.first) {
                return left.second < right.second;
            } else {
                return false;
            }
        }
    }; // struct Less
}; // struct IndexPairTraits


} // namespace detail


// @todo change to a more efficient set
template <class TIndex>
using NeighborSet = std::set<TIndex>;


/// @brief Collect all edges of an undirected graph.
template <class TIndex, class TValue>
std::vector<NeighborSet<TIndex>> collectNeighbors(const CSRMatrix<TIndex,TValue>& rMatrix)
{
    std::vector<NeighborSet<TIndex>> neighbors(rMatrix.columnCount);

    #pragma omp parallel
    {
        #pragma omp for
        for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
            const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
            const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

            for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
                const TValue value = rMatrix.pNonzeros[iEntry];
                if (value && iRow != iColumn) {
                    neighbors[iRow].insert(iColumn);
                }
            } // for iEntry in range(iRowBegin, iRowEnd)
        } // for iRow in range(rowCount)
    } // omp parallel

    return neighbors;
}


template <class TIndex, class TValue, class TColor>
int Color(const CSRMatrix<TIndex,TValue>& rMatrix,
          TColor* pColors,
          const ColoringSettings settings)
{
    // Cheap sanity checks
    if (rMatrix.rowCount < 0)                    return MCGS_FAILURE;
    if (rMatrix.nonzeroCount < 0)                return MCGS_FAILURE;
    if (!rMatrix.pRowExtents)                    return MCGS_FAILURE;
    if (!rMatrix.pColumnIndices)                 return MCGS_FAILURE;
    if (!rMatrix.pNonzeros)                      return MCGS_FAILURE;
    if (!pColors)                                return MCGS_FAILURE;
    if (rMatrix.rowCount != rMatrix.columnCount) return MCGS_FAILURE;

    // Collect all edges of the graph
    // (symmetric version of the input matrix)
    std::vector<NeighborSet<TIndex>> neighbors = collectNeighbors(rMatrix);

    // Find the minimum and maximum vertex degrees.
    TIndex minDegree = std::numeric_limits<TIndex>::max();
    TIndex maxDegree = 0;

    #pragma omp parallel for reduction(min: minDegree) reduction(max: maxDegree)
    for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        const TIndex degree = static_cast<TIndex>(neighbors[iRow].size());
        minDegree = std::min(minDegree, degree);
        maxDegree = std::max(maxDegree, degree);
    } // for iRow in range(rowCount)

    std::cout << "max vertex degree: " << maxDegree << "\nmin vertex degree: " << minDegree << std::endl;

    // Allocate the palette of every vertex to the max possible (maximum vertex degree).
    // An extra entry at the end of each palette indicates the palette's actual size.
    std::vector<TColor> palettes(rMatrix.columnCount * (maxDegree + 1));
    const TIndex shrunkPaletteSize = std::max(
        TIndex(1),
        TIndex(double(maxDegree) / double(std::max(TIndex(1), minDegree))));

    // Initialize the palette of all vertices to a shrunk set
    #pragma omp parallel for
    for (TIndex iVertex=0; iVertex<rMatrix.columnCount; ++iVertex) {
        const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
        const auto itPaletteEnd = itPaletteBegin + shrunkPaletteSize;
        std::iota(itPaletteBegin, itPaletteEnd, TColor(0));
        *(itPaletteBegin + maxDegree) = static_cast<TColor>(shrunkPaletteSize);
    } // for itPaletteBegin

    // Track vertices that need to be colored.
    std::vector<TIndex> toVisit;
    toVisit.reserve(rMatrix.columnCount);
    for (TIndex iVertex=0ul; iVertex<rMatrix.columnCount; ++iVertex) {
        toVisit.push_back(iVertex);
    }

    // Keep coloring until all vertices are colored.
    std::mt19937 randomGenerator(0);
    using UniformDistribution = std::uniform_int_distribution<TColor>;
    std::size_t iterationCount = 0ul;

    while (!toVisit.empty()) {
        std::cout << "coloring iteration " << iterationCount++ << " (" << toVisit.size() << " left to color)" << std::endl;
        // Assign random colors to each remaining vertex from their palette.
        #pragma omp parallel for firstprivate(randomGenerator)
        for (std::size_t iVisit=0ul; iVisit<toVisit.size(); ++iVisit) {
            const TIndex iVertex = toVisit[iVisit];
            const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
            const TColor paletteSize = *(itPaletteBegin + maxDegree);
            const TColor colorIndex = UniformDistribution(TColor(0), paletteSize)(randomGenerator);
            pColors[iVertex] = itPaletteBegin[colorIndex];
        }

        // Check for conflicts and remove colored vertices.
        std::vector<TIndex> toErase;
        for (TIndex iVertex : toVisit) {
            const TColor currentColor = pColors[iVertex];
            const auto neighborCount = static_cast<TIndex>(neighbors[iVertex].size());

            // If there's only one conflict, keep the coloring of the vertex with the higher index.
            bool conflict = false;
            bool colored = true;

            for (TIndex iNeighbor : neighbors[iVertex]) {
                const TColor neighborColor = pColors[iNeighbor];
                if (neighborColor == currentColor) {
                    if (conflict) {
                        // Multiple conflicts => give up on this vertex.
                        colored = false;
                        break;
                    } else if (iVertex < iNeighbor) {
                        // This is the first conflict, but the current vertex
                        // has a lower index than the neighbor it's in conflict with
                        // => give up on this vertex.
                        conflict = true;
                        colored = false;
                        break;
                    } else {
                        // This is the first conflict and the current vertex
                        // has a higher index than its conflicting neighbor,
                        // winning the tiebreaker => hang on to this vertex.
                        conflict = true;
                    }
                }
            } // for itRow in range(itRowBegin, itRowEnd)

            // If the current vertex has a valid color, remove it
            // from the remaining set and forbid its neighbors
            // and remove its color from the palettes of its neighbors.
            if (colored) {
                toErase.push_back(iVertex);
            }
        } // for iVertex in toVisit

        if (toErase.empty()) {
            // Failed to color any vertices => extend the palette of some random vertices
            const TIndex iVertex = UniformDistribution(0, toErase.size())(randomGenerator);
            const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
            TColor& rPaletteSize = *(itPaletteBegin + maxDegree);
            const TColor newColor = *(std::max_element(
                itPaletteBegin,
                itPaletteBegin + maxDegree
            )) + 1;
            itPaletteBegin[++rPaletteSize] = newColor;
        } else {
            // Remove colored vertices from the remaining set, as well
            // as their color from the palette of their neighbors.
            std::sort(toErase.begin(), toErase.end());
            {
                std::vector<TIndex> tmp;
                std::set_difference(toVisit.begin(), toVisit.end(),
                                    toErase.begin(), toErase.end(),
                                    std::back_inserter(tmp));
                toVisit.swap(tmp);
            }

            for (auto iVertex : toErase) {
                for (TIndex iNeighbor : neighbors[iVertex]) {
                    const auto itPaletteBegin = palettes.begin() + iNeighbor * (maxDegree + 1);
                    TColor& rPaletteSize = itPaletteBegin[maxDegree + 1];
                    const auto itPaletteEnd = itPaletteBegin + rPaletteSize;

                    if (std::remove(itPaletteBegin, itPaletteEnd, pColors[iVertex]) != itPaletteEnd) --rPaletteSize;

                    if (!rPaletteSize) {
                        // The neighbor's palette ran out of colors => find the
                        // highest color it ever had and add one higher to the
                        // empty palette.
                        const TColor newColor = *(std::max_element(
                            itPaletteBegin,
                            itPaletteBegin + maxDegree
                        )) + 1;
                        *itPaletteBegin = newColor;
                        ++rPaletteSize;
                    }
                } // for itRow in range(itRowBegin, itRowEnd)
            } // for iVertex in toErase
        }
    } // while toVisit

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_COLORING(TIndex, TValue, TColor)       \
    template int Color(const CSRMatrix<TIndex,TValue>& rMatrix, \
                       TColor* pColors,                         \
                       const ColoringSettings settings);

MCGS_INSTANTIATE_COLORING(int, double, unsigned);

MCGS_INSTANTIATE_COLORING(long, double, unsigned);

MCGS_INSTANTIATE_COLORING(unsigned, double, unsigned);

MCGS_INSTANTIATE_COLORING(std::size_t, double, unsigned);

MCGS_INSTANTIATE_COLORING(std::size_t, double, std::size_t);

#undef MCGS_INSTANTIATE_COLORING


} // namespace mcgs
