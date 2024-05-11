// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::solve, mcgs::CSRAdaptor

// --- STL Includes ---
#include <cstddef> // std::size_t


namespace mcgs {


template <class TIndex, class TValue, class TColor>
int solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const Partition<TIndex,TColor>* pPartition,
          const SolveSettings<TIndex,TValue> settings)
{
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
