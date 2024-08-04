// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::color
#include "multithreading.hpp"

// --- STL Includes ---
#include <vector> // std::vector
#include <cstddef> // std::size_t
#include <algorithm> // std::min, std::max, std::equal_range
#include <numeric> // std::iota
#include <random> // std::mt19937, std::uniform_int_distribution
#include <iostream> // std::cout, std::cerr


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


template <class TIndex>
using NeighborSet = std::vector<TIndex>;


/// @brief Collect all edges of an undirected graph.
template <class TIndex, class TValue>
std::vector<NeighborSet<TIndex>> collectNeighbors(const CSRAdaptor<TIndex,TValue>& rMatrix,
                                                  const ColorSettings<TIndex,TValue> settings,
                                                  [[maybe_unused]] MCGS_MUTEX_ARRAY& rMutexes)
{
    std::vector<NeighborSet<TIndex>> neighbors(rMatrix.columnCount);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
        const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            const TValue value = rMatrix.pNonzeros[iEntry];
            if (settings.tolerance < std::abs(value) && iRow != iColumn) {
                {
                    MCGS_ACQUIRE_MUTEX(rMutexes[iRow]);
                    const auto [itBegin, itEnd] = std::equal_range(neighbors[iRow].begin(),
                                                                   neighbors[iRow].end(),
                                                                   iColumn);
                    if (itBegin == itEnd) neighbors[iRow].insert(itBegin, iColumn);
                    MCGS_RELEASE_MUTEX(rMutexes[iRow]);
                }

                {
                    MCGS_ACQUIRE_MUTEX(rMutexes[iColumn]);
                    const auto [itBegin, itEnd] = std::equal_range(neighbors[iColumn].begin(),
                                                                   neighbors[iColumn].end(),
                                                                   iRow);
                    if (itBegin == itEnd) neighbors[iColumn].insert(itBegin, iRow);
                    MCGS_RELEASE_MUTEX(rMutexes[iColumn]);
                }
            }
        } // for iEntry in range(iRowBegin, iRowEnd)
    } // for iRow in range(rowCount)

    return neighbors;
}


using Mask = std::vector<char>;


template <class TIndex, class TColor>
bool isColored(const TIndex iVertex,
               const NeighborSet<TIndex>* pNeighborMap,
               const TColor* pColors,
               const Mask& rColoredMask)
{
    const TColor currentColor = pColors[iVertex];

    // If there's only one conflict, keep the coloring of the vertex with the higher index.
    bool colored = true;

    for (const TIndex iNeighbor : pNeighborMap[iVertex]) {
        const TColor neighborColor = pColors[iNeighbor];
        if (neighborColor == currentColor) {
            if (iVertex < iNeighbor) {
                // The current vertex has a lower index than the neighbor
                // it's in conflict with => give up on this vertex.
                colored = false;
                break;
            } else if (rColoredMask[iNeighbor]) {
                // Although the current vertex would win a tiebreaker against
                // its neighbor it's in conflict with, the neighbor's color is
                // already set and cannot be changed => give up on this vertex.
                colored = false;
                break;
            }

            // Otherwise, the current vertex has a higher index
            // than its conflicting neighbor, who is still waiting
            // to be colored, winning the tiebreaker
            // => hang on to this vertex.
        }
    } // for iNeighbor in neighbors[iVertex]

    return colored;
}


template <class TColor>
struct Palette
{
    TColor maxColor;

    std::vector<TColor> palette;
}; // struct Palette


template <class TColor>
bool extendPalette(Palette<TColor>& rPalette)
{
    rPalette.palette.push_back(++rPalette.maxColor);
    return true;
}


template <class TColor>
void removeFromPalette(const TColor color,
                       Palette<TColor>& rPalette)
{
    // The palette's colors are assumed to be sorted
    const auto [itBegin, itEnd] = std::equal_range(rPalette.palette.begin(),
                                                   rPalette.palette.end(),
                                                   color);
    rPalette.palette.erase(itBegin, itEnd);
}


template <class TIndex, class TValue, class TColor>
int color(TColor* pColors,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const ColorSettings<TIndex,TValue> settings)
{
    // Cheap sanity checks
    if (rMatrix.rowCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of rows " << rMatrix.rowCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.columnCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of columns " << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.nonzeroCount < 0) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: invalid number of nonzeros " << rMatrix.nonzeroCount << "\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pRowExtents) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing row data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pColumnIndices) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing column data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pNonzeros) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing nonzeros\n";
        return MCGS_FAILURE;
    }

    if (!pColors) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: missing output array\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.rowCount != rMatrix.columnCount) {
        if (1 < settings.verbosity) std::cerr << "mcgs: error: expecting a square matrix, but got "
                                              << rMatrix.rowCount << "x" << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

    MCGS_MUTEX_ARRAY mutexes(rMatrix.rowCount);
    for ([[maybe_unused]] MCGS_MUTEX& rMutex : mutexes) MCGS_INITIALIZE_MUTEX(rMutex);

    // Collect all edges of the graph
    // (symmetric version of the input matrix)
    const auto neighbors = collectNeighbors(rMatrix, settings, mutexes);

    // Find the minimum and maximum vertex degrees.
    TIndex minDegree = std::numeric_limits<TIndex>::max();
    TIndex maxDegree = 0;

    #ifdef MCGS_OPENMP
    #pragma omp parallel for reduction(min: minDegree) reduction(max: maxDegree)
    #endif
    for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        const TIndex degree = static_cast<TIndex>(neighbors[iRow].size());
        minDegree = std::min(minDegree, degree);
        maxDegree = std::max(maxDegree, degree);
    } // for iRow in range(rowCount)

    if (2 <= settings.verbosity) {
        std::cout << "mcgs: max vertex degree is " << maxDegree << '\n'
                  << "mcgs: min vertex degree is " << minDegree << '\n'
                  ;
    }

    // Allocate the palette of every vertex to the max possible (maximum vertex degree).
    // An extra entry at the end of each palette indicates the palette's actual size.
    const TIndex shrinkingFactor = 0 < settings.shrinkingFactor ?
                                   static_cast<TIndex>(settings.shrinkingFactor) :
                                   std::max(TIndex(1), minDegree);

    const TIndex initialPaletteSize = std::max(
        TIndex(1),
        TIndex(double(maxDegree) / double(shrinkingFactor)));

    if (3 <= settings.verbosity) {
        std::cout << "mcgs: initial palette size is " << initialPaletteSize << std::endl;
    }

    // Initialize the palette of all vertices to a shrunk set
    std::vector<Palette<TColor>> palettes(rMatrix.columnCount);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (TIndex iVertex=0; iVertex<rMatrix.columnCount; ++iVertex) {
        palettes[iVertex].palette.resize(initialPaletteSize);
        std::iota(palettes[iVertex].palette.begin(), palettes[iVertex].palette.end(), TColor(0));
        palettes[iVertex].maxColor = initialPaletteSize;
    } // for iVertex in range(columnCount)

    // Track vertices that need to be colored.
    Mask coloredMask(rMatrix.columnCount, false);
    std::vector<TIndex> uncolored;

    // Keep coloring until all vertices are colored.
    std::size_t iterationCount = 0ul;
    int stallCounter = 0;

    do {
        uncolored.clear();
        for (TIndex iVertex=0; iVertex<static_cast<TIndex>(coloredMask.size()); ++iVertex) {
            if (!coloredMask[iVertex]) uncolored.push_back(iVertex);
        }
        const TIndex uncoloredCount = uncolored.size();

        if (3 <= settings.verbosity) {
            std::cout << "mcgs: coloring iteration " << iterationCount++
                      << " (" << uncoloredCount << "/" << rMatrix.columnCount
                      << " left to color)\n";
        }

        // Assign random colors to each remaining vertex from their palette.
        #ifdef MCGS_OPENMP
        #pragma omp parallel
        #endif
        {

            #ifdef MCGS_OPENMP
                const int seed = omp_get_thread_num();
            #else
                const int seed = 0;
            #endif
            std::mt19937 randomGenerator(seed);



            #ifdef MCGS_OPENMP
            #pragma omp for
            #endif
            for (TIndex iVisit=0ul; iVisit<uncoloredCount; ++iVisit) {
                const TIndex iVertex = uncolored[iVisit];
                const TColor paletteSize = palettes[iVertex].palette.size();
                const TColor iColorMax = paletteSize ? paletteSize - 1 : static_cast<TColor>(0);
                const TColor iColor = std::uniform_int_distribution<TColor>(TColor(0), iColorMax)(randomGenerator);
                pColors[iVertex] = palettes[iVertex].palette[iColor];
            }

            #ifdef MCGS_OPENMP
            #pragma omp for
            #endif
            for (TIndex iVisit=0; iVisit<uncoloredCount; ++iVisit) {
                const TIndex iVertex = uncolored[iVisit];
                const bool colored = isColored(iVertex, neighbors.data(), pColors, coloredMask);

                // If the current vertex has a valid color, remove it
                // from the remaining set and remove its color from the
                // palettes of its neighbors.
                if (colored) {
                    for (TIndex iNeighbor : neighbors[iVertex]) {
                        if (!coloredMask[iNeighbor]) {
                            MCGS_ACQUIRE_MUTEX(mutexes[iNeighbor]);
                            removeFromPalette(pColors[iVertex], palettes[iNeighbor]);
                            if (palettes[iNeighbor].palette.empty() && palettes[iNeighbor].maxColor + 1 != pColors[iVertex]) {
                                extendPalette(palettes[iNeighbor]);
                            }
                            MCGS_RELEASE_MUTEX(mutexes[iNeighbor]);
                        } // if !coloredMask[iNeighbor]
                    } // for iNeighbor in neighbors[iVertex]

                    coloredMask[iVertex] = true;
                } // if colored
            } // for iVertex in uncolored
        } // omp parallel

        if (!uncolored.empty() && uncolored.size() == static_cast<std::size_t>(std::count_if(coloredMask.begin(), coloredMask.end(), [](const auto flag) {return !flag;}))) {
            // Failed to color any vertices => extend the palette of some random vertices
            const TIndex maxExtensions = std::max(TIndex(1), TIndex(5 * uncolored.size() / 100));
            TIndex extensionCounter = 0;

            #ifdef MCGS_OPENMP
            #pragma omp parallel for reduction(+: extensionCounter)
            #endif
            for (TIndex iUncolored=0; iUncolored<uncoloredCount; ++iUncolored) {
                if (extensionCounter < maxExtensions && extendPalette(palettes[uncolored[iUncolored]])) {
                    ++extensionCounter;
                }
            }

            if (!extensionCounter) {
                ++stallCounter;
                if (0 <= settings.maxStallCount && settings.maxStallCount <= stallCounter) {
                    if (1 <= settings.verbosity) std::cerr << "mcgs: error: reached stall limit (" << settings.maxStallCount << ")\n";
                    return MCGS_FAILURE;
                }
            } else {
                stallCounter = 0;
            }
        } else {
            stallCounter = 0;
        }
    } while (!uncolored.empty());

    for ([[maybe_unused]] MCGS_MUTEX& rMutex : mutexes) MCGS_DEINITIALIZE_MUTEX(rMutex);

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_COLOR(TIndex, TValue, TColor)              \
    template int color(TColor* pColors,                             \
                       const CSRAdaptor<TIndex,TValue>& rMatrix,    \
                       const ColorSettings<TIndex,TValue> settings);

MCGS_INSTANTIATE_COLOR(int, double, unsigned);
MCGS_INSTANTIATE_COLOR(long, double, unsigned);
MCGS_INSTANTIATE_COLOR(unsigned, double, unsigned);
MCGS_INSTANTIATE_COLOR(std::size_t, double, unsigned);
MCGS_INSTANTIATE_COLOR(std::size_t, double, std::size_t);

#undef MCGS_INSTANTIATE_COLOR


} // namespace mcgs
