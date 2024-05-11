// --- External Includes ---
#include <omp.h>

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::GaussSeidel

// --- STL Includes ---
#include <vector> // std::vector
#include <set> // std::set, std::pair
#include <unordered_set> // std::unordered_set
#include <unordered_map> // std::unordered_map
#include <cstddef> // std::size_t
#include <algorithm> // std::min, std::max
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


// @todo change to a more efficient set
template <class TIndex>
using NeighborSet = std::set<TIndex>;


/// @brief Collect all edges of an undirected graph.
template <class TIndex, class TValue>
std::vector<NeighborSet<TIndex>> collectNeighbors(const CSRMatrix<TIndex,TValue>& rMatrix)
{
    std::vector<NeighborSet<TIndex>> neighbors(rMatrix.columnCount);
    std::vector<omp_lock_t> locks;
    locks.reserve(neighbors.size());
    for (std::size_t iLock=0ul; iLock<neighbors.size(); ++iLock) {
        locks.emplace_back();
        omp_init_lock(&locks.back());
    }

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
                    omp_set_lock(&locks[std::min(iRow, iColumn)]);
                    omp_set_lock(&locks[std::max(iRow, iColumn)]);
                    neighbors[iRow].insert(iColumn);
                    neighbors[iColumn].insert(iRow);
                    omp_unset_lock(&locks[std::max(iRow, iColumn)]);
                    omp_unset_lock(&locks[std::min(iRow, iColumn)]);
                }
            } // for iEntry in range(iRowBegin, iRowEnd)
        } // for iRow in range(rowCount)
    } // omp parallel

    for (auto& rLock : locks) {
        omp_destroy_lock(&rLock);
    }

    return neighbors;
}


template <class TIndex, class TValue, class TColor>
int Color(const CSRMatrix<TIndex,TValue>& rMatrix,
          TColor* pColors,
          const ColoringSettings settings)
{
    // Cheap sanity checks
    if (rMatrix.rowCount < 0) {
        if (1 < settings.verbosity) std::cerr << "Error: invalid number of rows " << rMatrix.rowCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.columnCount < 0) {
        if (1 < settings.verbosity) std::cerr << "Error: invalid number of columns " << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.nonzeroCount < 0) {
        if (1 < settings.verbosity) std::cerr << "Error: invalid number of nonzeros " << rMatrix.nonzeroCount << "\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pRowExtents) {
        if (1 < settings.verbosity) std::cerr << "Error: missing row data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pColumnIndices) {
        if (1 < settings.verbosity) std::cerr << "Error: missing column data\n";
        return MCGS_FAILURE;
    }

    if (!rMatrix.pNonzeros) {
        if (1 < settings.verbosity) std::cerr << "Error: missing nonzeros\n";
        return MCGS_FAILURE;
    }

    if (!pColors) {
        if (1 < settings.verbosity) std::cerr << "Error: missing output array\n";
        return MCGS_FAILURE;
    }

    if (rMatrix.rowCount != rMatrix.columnCount) {
        if (1 < settings.verbosity) std::cerr << "Error: expecting a square matrix, but got "
                                              << rMatrix.rowCount << "x" << rMatrix.columnCount << "\n";
        return MCGS_FAILURE;
    }

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

    if (2 <= settings.verbosity) {
        std::cout << "max vertex degree: " << maxDegree
                  << "\nmin vertex degree: " << minDegree
                  << "\n";
    }

    // Allocate the palette of every vertex to the max possible (maximum vertex degree).
    // An extra entry at the end of each palette indicates the palette's actual size.
    std::vector<TColor> palettes(rMatrix.columnCount * (maxDegree + 1));
    const TIndex shrinkingFactor = 0 < settings.shrinkingFactor ?
                                   static_cast<TIndex>(settings.shrinkingFactor) :
                                   std::max(TIndex(1), minDegree);

    const TIndex initialPaletteSize = std::max(
        TIndex(1),
        TIndex(double(maxDegree) / double(shrinkingFactor)));

    if (3 <= settings.verbosity) {
        std::cout << "initial palette size: " << initialPaletteSize << std::endl;
    }

    // Initialize the palette of all vertices to a shrunk set
    #pragma omp parallel for
    for (TIndex iVertex=0; iVertex<rMatrix.columnCount; ++iVertex) {
        const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
        const auto itPaletteEnd = itPaletteBegin + initialPaletteSize;
        std::iota(itPaletteBegin, itPaletteEnd, TColor(0));
        std::fill(itPaletteEnd, itPaletteBegin + maxDegree, TColor(0));
        itPaletteBegin[maxDegree] = static_cast<TColor>(initialPaletteSize);
    } // for itPaletteBegin

    // Track vertices that need to be colored.
    std::vector<TIndex> toVisit;
    toVisit.reserve(rMatrix.columnCount);
    for (TIndex iVertex=0ul; iVertex<rMatrix.columnCount; ++iVertex) {
        toVisit.push_back(iVertex);
    }

    // Keep coloring until all vertices are colored.
    using UniformDistribution = std::uniform_int_distribution<TColor>;
    std::size_t iterationCount = 0ul;

    while (!toVisit.empty()) {
        const std::size_t visitCount = toVisit.size();

        if (3 <= settings.verbosity) {
            std::cout << "coloring iteration " << iterationCount++
                      << " (" << visitCount << "/" << rMatrix.columnCount
                      << " left to color)\n";
        }

        // Assign random colors to each remaining vertex from their palette.
        #pragma omp parallel
        {
            std::mt19937 randomGenerator(omp_get_thread_num());

            #pragma omp for
            for (std::size_t iVisit=0ul; iVisit<toVisit.size(); ++iVisit) {
                const TIndex iVertex = toVisit[iVisit];
                const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
                const TColor paletteSize = itPaletteBegin[maxDegree];
                const TColor colorIndex = UniformDistribution(TColor(0), paletteSize)(randomGenerator);
                pColors[iVertex] = itPaletteBegin[colorIndex];
            }
        } // omp parallel

        // Check for conflicts and remove colored vertices.
        std::unordered_map<TIndex,omp_lock_t> locks;
        for (auto iVertex : toVisit) {
            omp_init_lock(&locks.emplace(iVertex, omp_lock_t()).first->second);
        }

        #pragma omp parallel
        {
            std::vector<TIndex> verticesToErase;

            #pragma omp for
            for (std::size_t iVisit=0; iVisit<toVisit.size(); ++iVisit) {
                const TIndex iVertex = toVisit[iVisit];
                const TColor currentColor = pColors[iVertex];

                // If there's only one conflict, keep the coloring of the vertex with the higher index.
                bool conflict = false;
                bool colored = true;

                for (const TIndex iNeighbor : neighbors[iVertex]) {
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
                            const auto [itEqualBegin, itEqualEnd] = std::equal_range(toVisit.begin(), toVisit.end(), iNeighbor);
                            if (itEqualBegin == itEqualEnd) {
                                // Although the current vertex would win a tiebreaker against
                                // its neighbor it's in conflict with, the neighbor's color is
                                // already set and cannot be changed => give up on this vertex.
                                conflict = true;
                                colored = false;
                                break;
                            } else {
                                // This is the first conflict and the current vertex
                                // has a higher index than its conflicting neighbor,
                                // who is still waiting to be colored, winning the tiebreaker
                                // => hang on to this vertex.
                                conflict = true;
                            }
                        }
                    }
                } // for iNeighbor in neighbors[iVertex]

                // If the current vertex has a valid color, remove it
                // from the remaining set and remove its color from the
                // palettes of its neighbors.
                if (colored) {
                    verticesToErase.push_back(iVertex);
                }
            } // for iVertex in toVisit

            //#pragma omp critical
            //{
            //    if (verticesToErase.size()) {
            //        std::cout << "erase "; for (auto i : verticesToErase) {std::cout << i << " ";}
            //        std::cout << std::endl;
            //    }
            //}

            #pragma omp for
            for (std::size_t iEntry=0; iEntry<verticesToErase.size(); ++iEntry) {
                const auto iVertex = verticesToErase[iEntry];
                for (TIndex iNeighbor : neighbors[iVertex]) {
                    const auto itNeighbor = locks.find(iNeighbor);

                    if (itNeighbor != locks.end()) {
                        omp_set_lock(&itNeighbor->second);

                        const auto itPaletteBegin = palettes.begin() + iNeighbor * (maxDegree + 1);
                        TColor& rPaletteSize = itPaletteBegin[maxDegree];
                        const auto itPaletteEnd = itPaletteBegin + rPaletteSize;

                        if (std::remove(itPaletteBegin, itPaletteEnd, pColors[iVertex]) != itPaletteEnd) {
                            --rPaletteSize;
                        }

                        if (!rPaletteSize) {
                            // The neighbor's palette ran out of colors => find the
                            // highest color it ever had and add one higher to the
                            // empty palette.
                            const TColor newColor = (*std::max_element(
                                itPaletteBegin,
                                itPaletteBegin + maxDegree
                            )) + 1;
                            *itPaletteBegin = newColor;
                            rPaletteSize = 1;
                        }

                        omp_unset_lock(&itNeighbor->second);
                    } // if itNeighbor
                } // for iNeighbor in neighbors[iVertex]
            } // for iVertex in toErase

            // Remove colored vertices from the remaining set, as well
            // as their color from the palette of their neighbors.
            if (!verticesToErase.empty()) {
                #pragma omp critical
                {
                    std::vector<TIndex> tmp;
                    std::set_difference(toVisit.begin(), toVisit.end(),
                                        verticesToErase.begin(), verticesToErase.end(),
                                        std::back_inserter(tmp));
                    toVisit.swap(tmp);
                }
            }

        } // omp parallel

        for (auto& [iVertex, rLock] : locks) {
            omp_destroy_lock(&rLock);
        }

        if (toVisit.size() == visitCount) {
            // Failed to color any vertices => extend the palette of some random vertices
            bool success = false;

            for (TIndex iVertex : toVisit) {
                const auto itPaletteBegin = palettes.begin() + iVertex * (maxDegree + 1);
                TColor& rPaletteSize = itPaletteBegin[maxDegree];

                if (rPaletteSize < maxDegree) {
                    // Find the largest color the vertex ever had, and add
                    // a color one greater to its palette.
                    const TColor newColor = (*std::max_element(
                        itPaletteBegin,
                        itPaletteBegin + maxDegree
                    )) + 1;
                    itPaletteBegin[rPaletteSize++] = newColor;
                    success = true;
                    break;
                }
            }

            if (!success) {
                if (1 <= settings.verbosity) std::cerr << "Error: all remaining nodes' palettes are full\n";
                return MCGS_FAILURE;
            }
        } // if toVisit.size() == visitCount
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