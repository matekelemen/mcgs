// --- Internal Includes ---
#include "mcgs/mcgs.hpp"
#include "parseMatrixMarket.hpp"

// --- STL Includes ---
#include <filesystem> // std::filesystem::path
#include <iostream> // std::cout, std::cerr
#include <fstream> // std::ifstream
#include <limits>


namespace mcgs {





} // namespace mcgs


int main(int argc, const char* const * argv)
{
    for (unsigned i_arg=1; i_arg<argc; ++i_arg) {
        std::filesystem::path matrixPath = argv[i_arg];
        mcgs::TestCSRMatrix matrix;

        {
            std::ifstream matrixFile(matrixPath);
            matrix = mcgs::parseMatrixMarket(matrixFile);
        }

        mcgs::CSRMatrix<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value> adaptor;
        adaptor.rowCount = matrix.rowCount;
        adaptor.columnCount = matrix.columnCount;
        adaptor.nonzeroCount = matrix.nonzeroCount;
        adaptor.pRowExtents = matrix.rowExtents.data();
        adaptor.pColumnIndices = matrix.columnIndices.data();
        adaptor.pNonzeros = matrix.nonzeros.data();
        std::cout << matrix.rowCount << "x" << matrix.columnCount << " matrix\n"
                  << "\t" << matrix.rowExtents.size()       << "\n"
                  << "\t" << matrix.columnIndices.size()    << "\n"
                  << "\t" << matrix.nonzeros.size()         << "\n";
        std::cout << std::endl;

        std::vector<unsigned> colors(matrix.columnCount, std::numeric_limits<unsigned>::max());
        mcgs::Color(adaptor, colors.data(), {});
    }

    return 0;
}
