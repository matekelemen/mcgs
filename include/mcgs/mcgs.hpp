#pragma once


#define MCGS_EXPORT_SYMBOL __attribute__((visibility ("default")))


namespace mcgs {


#define MCGS_SUCCESS 0
#define MCGS_FAILURE 1


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


template <class TIndex, class TValue>
struct ColorSettings
{
    int shrinkingFactor;
    int maxStallCount;
    TValue tolerance;
    unsigned short verbosity;

    ColorSettings() noexcept
        : shrinkingFactor(-1),
          maxStallCount(3),
          tolerance(0),
          verbosity(1)
    {}
};


template <class TIndex, class TValue, class TColor>
MCGS_EXPORT_SYMBOL
int color(TColor* pColors,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const ColorSettings<TIndex,TValue> settings);


template <class TIndex>
class Partition;


template <class TIndex, class TColor>
MCGS_EXPORT_SYMBOL
[[nodiscard]] Partition<TIndex>* makePartition(const TColor* pColors,
                                               const TIndex columnCount);


template <class TIndex>
MCGS_EXPORT_SYMBOL
void destroyPartition(Partition<TIndex>* pPartition);


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
[[nodiscard]] Partition<TIndex>* reorder(const TIndex rowCount, const TIndex columnCount, const TIndex nonzeroCount,
                                         TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
                                         TValue* pRHS,
                                         const Partition<TIndex>* pPartition);


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int revertReorder(TValue* pRHS, const Partition<TIndex>* pPartition);


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
TValue residual(const CSRAdaptor<TIndex,TValue>& rMatrix,
                const TValue* pSolution,
                const TValue* pRHS) noexcept;


enum struct Parallelization
{
    None        = 1,
    RowWise     = 2,
    EntryWise   = 4
};


template <class TIndex, class TValue>
struct SolveSettings
{
    TValue residualAbsoluteTolerance;
    TValue residualRelativeTolerance;
    TIndex maxIterations;
    TValue relaxation;
    Parallelization parallelization;
    unsigned short verbosity;

    SolveSettings() noexcept
        : residualAbsoluteTolerance(-1),
          residualRelativeTolerance(-1),
          maxIterations(1),
          relaxation(1),
          parallelization(Parallelization::RowWise),
          verbosity(1)
    {}
};


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const SolveSettings<TIndex,TValue> settings);


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex>* pPartition,
          const SolveSettings<TIndex,TValue> settings);


} // namespace mcgs


#undef MCGS_EXPORT_SYMBOL
