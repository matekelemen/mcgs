#pragma once

namespace mcgs {


#define MCGS_SUCCESS 0
#define MCGS_FAILURE 1


struct ColoringSettings
{
    int shrinkingFactor;

    unsigned short verbosity;

    ColoringSettings() noexcept
        : shrinkingFactor(-1),
          verbosity(0)
    {}
};


template <class TIndex, class TValue>
struct SolutionSettings
{
    TValue residualAbsoluteTolerance;
    TValue residualRelativeTolerance;
    TIndex maxIterations;

    SolutionSettings() noexcept
        : residualAbsoluteTolerance(1e-6),
          residualRelativeTolerance(1e-4),
          maxIterations(1e3)
    {}
};


template <class TIndex, class TValue>
struct CSRMatrix
{
    TIndex rowCount;
    TIndex columnCount;
    TIndex nonzeroCount;
    TIndex* pRowExtents;
    TIndex* pColumnIndices;
    TValue* pNonzeros;
};


template <class TIndex, class TValue, class TColor>
int Color(const CSRMatrix<TIndex,TValue>& rMatrix,
          TColor* pColors,
          const ColoringSettings settings);


template <class TIndex, class TValue, class TColor>
int Solve(const CSRMatrix<TIndex,TValue>& rMatrix,
          const TColor* pColors,
          TValue* pSolution,
          const SolutionSettings<TIndex,TValue> settings);


} // namespace mcgs
