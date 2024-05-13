// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::solve, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition

// --- STL Includes ---
#include <cstddef> // std::size_t
#include <omp.h>
#include <vector> // std::vector
#include <algorithm> // std::copy
#include <cmath> // std::sqrt
#include <iostream> // std::cout, std::cerr


namespace mcgs {


bool isWorthParallelizing(const std::size_t independentRowCount,
                          const std::size_t threadCount)
{
    return threadCount < 12 * independentRowCount;
}


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


template <class TIndex, class TValue, class TColor>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex,TColor>* pPartition,
          const SolveSettings<TIndex,TValue> settings)
{
    // Serial version
    //return solve(pSolution, rMatrix, pRHS, settings);
    const auto threadCount = omp_get_num_threads();

    std::vector<TValue> buffer(rMatrix.columnCount);
    const TValue initialResidual = residual(rMatrix, pSolution, pRHS, buffer.data());

    for (TIndex iIteration=0; iIteration<settings.maxIterations; ++iIteration) {
        for (TIndex iPartition=0; iPartition<pPartition->size(); ++iPartition) {
            const auto itPartitionBegin = pPartition->begin(iPartition);
            const auto partitionSize = pPartition->size(iPartition);

            #define MCGS_SWEEP                                                                  \
                for (TIndex iLocal=0; iLocal<partitionSize; ++iLocal) {                         \
                    const TIndex iRow = itPartitionBegin[iLocal];                               \
                    TValue value = pRHS[iRow];                                                  \
                    TValue diagonal = 1;                                                        \
                                                                                                \
                    const TIndex iRowBegin = rMatrix.pRowExtents[iRow];                         \
                    const TIndex iRowEnd = rMatrix.pRowExtents[iRow + 1];                       \
                                                                                                \
                    for (TIndex iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {                   \
                        const TIndex iColumn = rMatrix.pColumnIndices[iEntry];                  \
                        const TValue product = rMatrix.pNonzeros[iEntry] * pSolution[iColumn];  \
                        if (iColumn == iRow) diagonal = rMatrix.pNonzeros[iEntry];              \
                        else value -= product;                                                  \
                    } /*for iEntry in range(iRowBegin, iRowEnd)*/                               \
                                                                                                \
                pSolution[iRow] += settings.relaxation * (value / diagonal - pSolution[iRow]);  \
            } // for iLocal in range()

            if (isWorthParallelizing(pPartition->size(iPartition), threadCount)) {
                #pragma omp parallel for
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


#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue, TColor)                          \
    template int solve<TIndex,TValue,TColor>(TValue*,                           \
                                             const CSRAdaptor<TIndex,TValue>&,  \
                                             const TValue*,                     \
                                             const Partition<TIndex,TColor>*,   \
                                             const SolveSettings<TIndex,TValue>);

MCGS_INSTANTIATE_SOLVE(int, double, unsigned);

MCGS_INSTANTIATE_SOLVE(long, double, unsigned);

MCGS_INSTANTIATE_SOLVE(unsigned, double, unsigned);

MCGS_INSTANTIATE_SOLVE(std::size_t, double, unsigned);

MCGS_INSTANTIATE_SOLVE(std::size_t, double, std::size_t);

#undef MCGS_INSTANTIATE_SOLVE


} // namespace mcgs
