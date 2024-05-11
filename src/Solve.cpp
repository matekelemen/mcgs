// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::Solve, mcgs::CSRAdaptor

// --- STL Includes ---
#include <cstddef> // std::size_t


namespace mcgs {


template <class TIndex, class TValue, class TColor>
int Solve(TValue* pSolution,
          const CSRAdaptor<TIndex,TValue>& rMatrix,
          const TValue* pRHS,
          const TColor* pColors,
          const SolveSettings<TIndex,TValue> settings)
{
    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue, TColor)                                          \
    template int Solve<TIndex,TValue,TColor>(TValue* pSolution,                                 \
                                              const CSRAdaptor<TIndex,TValue>& rMatrix,         \
                                              const TValue* pRHS,                               \
                                              const TColor* pColors,                            \
                                              const SolveSettings<TIndex,TValue> settings);

MCGS_INSTANTIATE_SOLVE(int, double, unsigned);

MCGS_INSTANTIATE_SOLVE(long, double, unsigned);

MCGS_INSTANTIATE_SOLVE(unsigned, double, unsigned);

MCGS_INSTANTIATE_SOLVE(std::size_t, double, unsigned);

MCGS_INSTANTIATE_SOLVE(std::size_t, double, std::size_t);

#undef MCGS_INSTANTIATE_SOLVE


} // namespace mcgs
