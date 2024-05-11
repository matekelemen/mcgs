#pragma once

namespace mcgs {


#define MCGS_SUCCESS 0
#define MCGS_FAILURE 1


struct ColorSettings
{
    int shrinkingFactor;

    unsigned short verbosity;

    ColorSettings() noexcept
        : shrinkingFactor(-1),
          verbosity(1)
    {}
};


template <class TIndex, class TValue>
struct SolveSettings
{
    TValue residualAbsoluteTolerance;
    TValue residualRelativeTolerance;
    TIndex maxIterations;

    SolveSettings() noexcept
        : residualAbsoluteTolerance(1e-6),
          residualRelativeTolerance(1e-4),
          maxIterations(1e3)
    {}
};


template <class TIndex, class TValue>
struct CSRAdaptor
{
    TIndex rowCount;
    TIndex columnCount;
    TIndex nonzeroCount;
    const TIndex* pRowExtents;
    const TIndex* pColumnIndices;
    const TValue* pNonzeros;

    CSRAdaptor() noexcept
        : rowCount(TIndex(0)),
          columnCount(TIndex(0)),
          nonzeroCount(TIndex(0)),
          pRowExtents(nullptr),
          pColumnIndices(nullptr),
          pNonzeros(nullptr)
    {}
};


template <class TIndex, class TColor>
struct Partition;


template <class TIndex, class TValue, class TColor>
int color(TColor* pColors,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const ColorSettings settings);


template <class TIndex, class TColor>
[[nodiscard]] Partition<TIndex,TColor>* makePartition(const TColor* pColors,
                                                      const TIndex columnCount);


template <class TIndex, class TColor>
void destroyPartition(Partition<TIndex,TColor>* pPartition);


template <class TIndex, class TValue, class TColor>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex,TColor>* pPartition,
          const SolveSettings<TIndex,TValue> settings);


} // namespace mcgs
