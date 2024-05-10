#include "mcgs/mcgs.hpp" // mcgs::Solve, mcgs::CSRMatrix
#include <cstddef> // std::size_t


namespace mcgs {


template <class TIndex, class TValue, class TColor>
void Solve(const CSRMatrix<TIndex,TValue>& rMatrix,
           const TColor* pColors,
           TValue* pSolution,
           const SolutionSettings<TIndex,TValue> settings)
{
}


#define MCGS_INSTANTIATE_SOLVE(TIndex, TValue, TColor)                                          \
    template void Solve<TIndex,TValue,TColor>(const CSRMatrix<TIndex,TValue>& rMatrix,          \
                                              const TColor* pColors,                            \
                                              TValue* pSolution,                                \
                                              const SolutionSettings<TIndex,TValue> settings);

MCGS_INSTANTIATE_SOLVE(int, double, unsigned);

MCGS_INSTANTIATE_SOLVE(long, double, unsigned);

MCGS_INSTANTIATE_SOLVE(unsigned, double, unsigned);

MCGS_INSTANTIATE_SOLVE(std::size_t, double, unsigned);

#undef MCGS_INSTANTIATE_SOLVE


} // namespace mcgs
