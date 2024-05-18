// --- External Includes ---
#ifdef MCGS_OPENMP
#include <omp.h> // omp_get_num_threads
#endif

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::solve, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <cstddef> // std::size_t
#include <vector> // std::vector
#include <algorithm> // std::copy, std::clamp
#include <cmath> // std::sqrt
#include <iostream> // std::cout, std::cerr


namespace mcgs {


template <class TIndex, class TValue>
TValue residual(const CSRAdaptor<TIndex,TValue>& rMatrix,
                const TValue* pSolution,
                const TValue* pRHS,
                TValue* buffer) noexcept
{
    std::copy(pRHS, pRHS + rMatrix.columnCount, buffer);
    TValue residual = 0;

    #ifdef MCGS_OPENMP
    #pragma omp parallel for reduction(+: residual)
    #endif
    for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        TValue& rResidualComponent = buffer[iRow];
        const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
        const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

        for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
            rResidualComponent -= rMatrix.pNonzeros[iEntry] * pSolution[iColumn];
        } // for iEntry in range(iRowBegin, iRowEnd)

        residual += rResidualComponent * rResidualComponent;
    } // for iRow in range(0, rowCount)

    return std::sqrt(residual);
}


template <class TIndex, class TValue>
int contiguousGaussSeidelSweep(TValue* pSolution,
                               const CSRAdaptor<TIndex,TValue>& rMatrix,
                               const TValue* pRHS,
                               const SolveSettings<TIndex,TValue> settings,
                               const TIndex iRowBegin,
                               const TIndex iRowEnd,
                               const int threadCount)
{
    if (iRowEnd < iRowBegin) {
        if (1 <= settings.verbosity) {
            std::cerr << "Error: invalid range [" << iRowBegin << ", " << iRowEnd << "[\n";
        }
        return MCGS_FAILURE;
    }

    #define MCGS_ROW_WISE_SWEEP                                                                 \
        for (TIndex iRow=iRowBegin; iRow<iRowEnd; ++iRow) {                                     \
            TValue value = pRHS[iRow];                                                          \
            TValue diagonal = 1;                                                                \
                                                                                                \
            const TIndex iEntryBegin = rMatrix.pRowExtents[iRow];                               \
            const TIndex iEntryEnd = rMatrix.pRowExtents[iRow + 1];                             \
                                                                                                \
            for (TIndex iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {                       \
                const TIndex iColumn = rMatrix.pColumnIndices[iEntry];                          \
                const TValue nonzero = rMatrix.pNonzeros[iEntry];                               \
                                                                                                \
                if (iColumn == iRow) diagonal = nonzero;                                        \
                else value -= nonzero * pSolution[iColumn];                                     \
            } /*for iEntry in range(iEntryBegin, iEntryEnd)*/                                   \
                                                                                                \
            pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);      \
        }

    if (1 < threadCount) {
        if (settings.parallelization == Parallelization::RowWise) {
            #ifdef MCGS_OPENMP
            #pragma omp parallel for
            #endif
            MCGS_ROW_WISE_SWEEP
        } /*if settings.parallelization == RowWise*/ else if (settings.parallelization == Parallelization::NonzeroWise) {
            const TIndex partitionRowCount = iRowEnd - iRowBegin;
            const auto itEntryBegin = rMatrix.pRowExtents + iRowBegin;
            const auto itEntryEnd = rMatrix.pRowExtents + iRowEnd;
            const TIndex iEntryBegin = *itEntryBegin;
            const TIndex iEntryEnd = *itEntryEnd;
            const TIndex entryCount = iEntryEnd - iEntryBegin;

            std::vector<TValue> diagonals(partitionRowCount);
            std::vector<TValue> updates(partitionRowCount);
            std::copy(pRHS + iRowBegin, pRHS + iRowEnd, updates.data());

            std::vector<TIndex> threadEntryExtents(threadCount + 1);
            threadEntryExtents.front() = iEntryBegin;
            {
                const TIndex chunkSize = entryCount / threadCount + (entryCount % threadCount ? 1 : 0);
                for (TIndex iEnd=1; iEnd<static_cast<TIndex>(threadCount) + 1; ++iEnd) {
                    threadEntryExtents[iEnd] = std::min(
                        iEntryEnd,
                        threadEntryExtents[iEnd - 1] + chunkSize
                    );
                }
            }

            #ifdef MCGS_OPENMP
            #pragma omp parallel
            #endif
            {
                std::vector<TValue> localUpdates(partitionRowCount, static_cast<TValue>(0));

                #ifdef MCGS_OPENMP
                const TIndex iThread = omp_get_thread_num();
                #else
                const TIndex iThread = 0;
                #endif

                const TIndex iLocalEntryBegin = threadEntryExtents[iThread];
                const TIndex iLocalEntryEnd = threadEntryExtents[iThread + 1];

                const auto itLocalEntryBegin = std::max(std::upper_bound(itEntryBegin, itEntryEnd, iLocalEntryBegin) - 1,
                                                        itEntryBegin);
                TIndex iRow = std::distance(rMatrix.pRowExtents, itLocalEntryBegin);
                TIndex iLocalRow = iRow - iRowBegin;
                const TIndex* itRowEnd = rMatrix.pRowExtents + iRow + 1;

                for (TIndex iEntry=iLocalEntryBegin; iEntry<iLocalEntryEnd; ++iEntry) {
                    while (*itRowEnd <= iEntry) {
                        ++iRow;
                        ++iLocalRow;
                        ++itRowEnd;
                    }

                    const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
                    const TValue nonzero = rMatrix.pNonzeros[iEntry];

                    if (iRow == iColumn) {
                        diagonals[iLocalRow] = nonzero;
                    } else {
                        localUpdates[iLocalRow] -= nonzero * pSolution[iColumn];
                    }
                } // for iEntry in range(iLocalEntryBegin, iLocalEntryEnd)

                #ifdef MCGS_OPENMP
                #pragma omp critical
                #endif
                {
                    for (TIndex iLocal=0; iLocal<static_cast<TIndex>(updates.size()); ++iLocal) {
                        updates[iLocal] += localUpdates[iLocal];
                    }
                } // omp critical
            } // omp parallel

            #ifdef MCGS_OPENMP
            #pragma omp parallel for
            #endif
            for (TIndex iRow=iRowBegin; iRow<iRowEnd; ++iRow) {
                const TIndex iLocalRow = iRow - iRowBegin;
                pSolution[iRow] += settings.relaxation * (updates[iLocalRow] / diagonals[iLocalRow] - pSolution[iRow]);
            }
        } /*if settings.parallelization == Parallelization::NonzeroWise*/ else if (settings.parallelization == Parallelization::None) {
            MCGS_ROW_WISE_SWEEP
        } // /*if settings.parallelization == Parallelization::None*/
    } /*if 1 < threadCount*/ else {
        MCGS_ROW_WISE_SWEEP
    }

    #undef MCGS_ROW_WISE_SWEEP
    return MCGS_SUCCESS;
}


template <class TIndex, class TValue>
int randomAccessGaussSeidelSweep(TValue* pSolution,
                                 const CSRAdaptor<TIndex,TValue>& rMatrix,
                                 const TValue* pRHS,
                                 const SolveSettings<TIndex,TValue> settings,
                                 const TIndex* itPartitionBegin,
                                 const TIndex* itPartitionEnd,
                                 const int threadCount)
{
    if (itPartitionEnd < itPartitionBegin) {
        if (1 <= settings.verbosity) {
            std::cerr << "Error: invalid partition range\n";
        }
        return MCGS_FAILURE;
    }

    if (settings.parallelization == Parallelization::NonzeroWise) {
        if (1 <= settings.verbosity) {
            std::cerr << "Error: cannot perform nonzero-wise parallelization on random access partitions\n";
        }
        return MCGS_FAILURE;
    }

    if (1 < threadCount && settings.parallelization != Parallelization::None) {
        const auto partitionSize = static_cast<typename Partition<TIndex>::size_type>(std::distance(itPartitionBegin, itPartitionEnd));

        #ifdef MCGS_OPENMP
        #pragma omp parallel for
        #endif
        for (typename Partition<TIndex>::size_type iLocal=0; iLocal<partitionSize; ++iLocal) {
            const TIndex iRow = itPartitionBegin[iLocal];
            TValue value = pRHS[iRow];
            TValue diagonal = 1;

            const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
            const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

            for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
                if (iColumn == iRow) diagonal = rMatrix.pNonzeros[iEntry];
                else value -= rMatrix.pNonzeros[iEntry] * pSolution[iColumn];
            } /*for iEntry in range(iRowBegin, iRowEnd)*/

            pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);
        } // for iLocal in range()
    } /*if 1 < threadCount && settings.parallelization != Parallelization::None*/ else {
        for (auto itRow=itPartitionBegin; itRow!=itPartitionEnd; ++itRow) {
            const TIndex iRow = *itRow;
            TValue value = pRHS[iRow];
            TValue diagonal = 1;

            const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
            const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

            for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
                if (iColumn == iRow) diagonal = rMatrix.pNonzeros[iEntry];
                else value -= rMatrix.pNonzeros[iEntry] * pSolution[iColumn];
            } /*for iEntry in range(iRowBegin, iRowEnd)*/

            pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);
        } // for iLocal in range()
    }
    return MCGS_SUCCESS;
}


template <class TIndex, class TValue>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const SolveSettings<TIndex,TValue> settings)
{
    if (settings.parallelization != Parallelization::None) {
        if (1 <= settings.verbosity) {
            std::cerr << "Error: parallel Gauss-Seidel requires a partition\n";
        }
        return MCGS_FAILURE;
    }

    std::vector<TValue> buffer(3 <= rMatrix.columnCount ? rMatrix.columnCount : 0);
    const TValue initialResidual = 3 <= settings.verbosity ?
                                   residual(rMatrix, pSolution, pRHS, buffer.data()) :
                                   static_cast<TValue>(1);

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        if (contiguousGaussSeidelSweep(pSolution,
                                       rMatrix,
                                       pRHS,
                                       settings,
                                       TIndex(0),
                                       rMatrix.rowCount,
                                       1) != MCGS_SUCCESS) {
            return MCGS_FAILURE;
        }

        if (3 <= settings.verbosity) {
            std::cout << "iteration " << iIteration
                      << " residual: "
                      << residual(rMatrix, pSolution, pRHS, buffer.data()) / initialResidual
                      << "\n";
        }
    }

    return MCGS_SUCCESS;
}


template <class TIndex, class TValue>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex>* pPartition,
          const SolveSettings<TIndex,TValue> settings)
{
    // Collect how many threads should execute each partition.
    #ifdef MCGS_OPENMP
    const int maxThreadCount = omp_get_max_threads();
    #else
    const int maxThreadCount = 1;
    #endif
    std::vector<int> threadCounts(pPartition->size());

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (int iPartition=0; iPartition<static_cast<int>(pPartition->size()); ++iPartition) {
        std::size_t nonzeroCount = 0ul;
        for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            const TIndex iRow = *itPartition;
            nonzeroCount += rMatrix.pRowExtents[iRow + 1] - rMatrix.pRowExtents[iRow];
        } // for itPartition in pPartition[iPartition]

        // @todo Find a dynamic way of approximating the optimal load of a single thread.
        //threadCounts[iPartition] = std::clamp(nonzeroCount / 1024,
        //                                      static_cast<std::size_t>(1),
        //                                      std::min(pPartition->size(iPartition),
        //                                               static_cast<std::size_t>(maxThreadCount)));
        threadCounts[iPartition] = maxThreadCount;
    }

    std::vector<TValue> buffer(rMatrix.columnCount);
    const TValue initialResidual = 3 <= settings.verbosity ?
                                   residual(rMatrix, pSolution, pRHS, buffer.data()) :
                                   static_cast<TValue>(1);

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
            const auto threadCount = threadCounts[iPartition];
            if (pPartition->isContiguous()) {
                if (contiguousGaussSeidelSweep(pSolution,
                                               rMatrix,
                                               pRHS,
                                               settings,
                                               *pPartition->begin(iPartition),
                                               *pPartition->end(iPartition),
                                               threadCount) != MCGS_SUCCESS) {
                    return MCGS_FAILURE;
                }
            } /*if pPartition->isContiguous()*/ else {
                if (randomAccessGaussSeidelSweep(pSolution,
                                                 rMatrix,
                                                 pRHS,
                                                 settings,
                                                 pPartition->begin(iPartition),
                                                 pPartition->end(iPartition),
                                                 threadCount) != MCGS_SUCCESS) {
                    return MCGS_FAILURE;
                }
            }
        } // for iPartition in range(partitionCount)

        if (3 <= settings.verbosity) {
            std::cout << "iteration " << iIteration
                      << " residual: "
                      << residual(rMatrix, pSolution, pRHS, buffer.data()) / initialResidual
                      << "\n";
        } // if 3 <= settings.verbosity
    } // for iIteration in range(settings.maxIterations)

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue)                              \
    template TValue residual(const CSRAdaptor<TIndex,TValue>&,              \
                             const TValue*,                                 \
                             const TValue*,                                 \
                             TValue*) noexcept;                             \
    template int solve<TIndex,TValue>(TValue*,                              \
                                      const CSRAdaptor<TIndex,TValue>&,     \
                                      const TValue*,                        \
                                      const SolveSettings<TIndex,TValue>);  \
    template int solve<TIndex,TValue>(TValue*,                              \
                                      const CSRAdaptor<TIndex,TValue>&,     \
                                      const TValue*,                        \
                                      const Partition<TIndex>*,             \
                                      const SolveSettings<TIndex,TValue>);

MCGS_INSTANTIATE_SOLVE(int, double);
MCGS_INSTANTIATE_SOLVE(long, double);
MCGS_INSTANTIATE_SOLVE(unsigned, double);
MCGS_INSTANTIATE_SOLVE(std::size_t, double);

#undef MCGS_INSTANTIATE_SOLVE


} // namespace mcgs
