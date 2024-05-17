// --- External Includes ---
#include <omp.h> // omp_get_num_threads

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

    #pragma omp parallel for reduction(+: residual)
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
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const SolveSettings<TIndex,TValue> settings)
{
    std::vector<TValue> buffer(rMatrix.columnCount);
    const TValue initialResidual = residual(rMatrix, pSolution, pRHS, buffer.data());

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        for (TIndex iRow=0; iRow<rMatrix.rowCount; ++iRow) {
            TValue value = pRHS[iRow];
            TValue diagonal = 1;

            const TIndex iRowBegin = rMatrix.pRowExtents[iRow];
            const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];

            for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const TIndex iColumn = rMatrix.pColumnIndices[iEntry];
                const TValue nonzero = rMatrix.pNonzeros[iEntry];

                if (iColumn == iRow) diagonal = nonzero;
                else value -= nonzero * pSolution[iColumn];
            } // for iEntry in range(iRowBegin, iRowEnd)

            pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);
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
    const auto maxThreadCount = omp_get_max_threads();
    std::vector<int> threadCounts(pPartition->size());

    #pragma omp parallel for
    for (int iPartition=0; iPartition<static_cast<int>(pPartition->size()); ++iPartition) {
        std::size_t nonzeroCount = 0ul;
        for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            const TIndex iRow = *itPartition;
            nonzeroCount += rMatrix.pRowExtents[iRow + 1] - rMatrix.pRowExtents[iRow];
        } // for itPartition in pPartition[iPartition]

        // @todo Find a dynamic way of approximating the optimal load of a single thread.
        threadCounts[iPartition] = std::clamp(nonzeroCount / 1024,
                                              static_cast<std::size_t>(1),
                                              std::min(pPartition->size(iPartition),
                                                       static_cast<std::size_t>(maxThreadCount)));
    }

    std::vector<TValue> buffer(rMatrix.columnCount);
    const TValue initialResidual = residual(rMatrix, pSolution, pRHS, buffer.data());

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
            const auto itPartitionBegin = pPartition->begin(iPartition);
            const auto partitionSize = pPartition->size(iPartition);

            #define MCGS_SWEEP                                                                          \
                for (typename Partition<TIndex>::size_type iLocal=0; iLocal<partitionSize; ++iLocal) {  \
                    const TIndex iRow = itPartitionBegin[iLocal];                                       \
                    TValue value = pRHS[iRow];                                                          \
                    TValue diagonal = 1;                                                                \
                                                                                                        \
                    const TIndex iRowBegin = rMatrix.pRowExtents[iRow];                                 \
                    const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];                               \
                                                                                                        \
                    for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {                           \
                        const TIndex iColumn = rMatrix.pColumnIndices[iEntry];                          \
                        if (iColumn == iRow) diagonal = rMatrix.pNonzeros[iEntry];                      \
                        else value -= rMatrix.pNonzeros[iEntry] * pSolution[iColumn];                   \
                    } /*for iEntry in range(iRowBegin, iRowEnd)*/                                       \
                                                                                                        \
                pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);          \
            } // for iLocal in range()

            const auto threadCount = threadCounts[iPartition];
            if (1 < threadCount) {
                #pragma omp parallel for num_threads(threadCount)
                MCGS_SWEEP
            } else {
                MCGS_SWEEP
            }

            #undef MCGS_SWEEP

        } // for iPartition in range(partitionCount)

        if (3 <= settings.verbosity) {
            std::cout << "iteration " << iIteration
                      << " residual: "
                      << residual(rMatrix, pSolution, pRHS, buffer.data()) / initialResidual
                      << "\n";
        }
    } // for iIteration in range(settings.maxIterations)

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue)                              \
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
