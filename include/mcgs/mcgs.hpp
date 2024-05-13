#pragma once

namespace mcgs {


#define MCGS_SUCCESS 0
#define MCGS_FAILURE 1


struct ColorSettings
{
    int shrinkingFactor;
    int maxStallCount;
    unsigned short verbosity;

    ColorSettings() noexcept
        : shrinkingFactor(-1),
          maxStallCount(3),
          verbosity(1)
    {}
};


template <class TIndex, class TValue>
struct SolveSettings
{
    TValue residualAbsoluteTolerance;
    TValue residualRelativeTolerance;
    TIndex maxIterations;
    TValue relaxation;
    unsigned short verbosity;

    SolveSettings() noexcept
        : residualAbsoluteTolerance(-1),
          residualRelativeTolerance(-1),
          maxIterations(1),
          relaxation(1),
          verbosity(1)
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
int reorder(const TIndex rowCount, const TIndex columnCount, const TIndex nonzeroCount,
            TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
            TValue* pRHS,
            Partition<TIndex,TColor>* pPartition);


template <class TIndex, class TValue, class TColor>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex,TColor>* pPartition,
          const SolveSettings<TIndex,TValue> settings);


} // namespace mcgs
