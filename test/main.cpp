// --- Internal Includes ---
#include "mcgs/mcgs.hpp"
#include "parseMatrixMarket.hpp"

// --- STL Includes ---
#include <filesystem> // std::filesystem::path
#include <iostream> // std::cout, std::cerr
#include <fstream> // std::ifstream
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <limits> // std::numeric_limits::max


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
        std::cout << std::endl;

        std::vector<unsigned> colors(matrix.columnCount, std::numeric_limits<unsigned>::max());
        mcgs::ColoringSettings settings;
        settings.verbosity = 1;
        settings.shrinkingFactor = 15;
        mcgs::Color(adaptor, colors.data(), settings);

        if (matrix.columnCount < 40) {
            std::cout << "colors: ";
            for (auto c : colors) std::cout << c << " ";
            std::cout << std::endl;
        }

        std::unordered_map<
            mcgs::TestCSRMatrix::Index,
            std::vector<mcgs::TestCSRMatrix::Index>
        > conflicts;
        std::unordered_set<unsigned> palette;
        for (std::size_t iRow=0; iRow<matrix.columnCount; ++iRow) {
            const auto currentColor = colors[iRow];
            palette.insert(currentColor);
            const auto iRowBegin = matrix.rowExtents[iRow];
            const auto iRowEnd = matrix.rowExtents[iRow + 1];
            for (std::size_t iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const mcgs::TestCSRMatrix::Index iColumn = matrix.columnIndices[iEntry];
                const auto neighborColor = colors[iColumn];
                if (matrix.nonzeros[iEntry] && iRow != iColumn) {
                    if (neighborColor == currentColor) {
                        conflicts.emplace(iRow, std::vector<mcgs::TestCSRMatrix::Index> {})
                            .first->second.push_back(iColumn);
                    }
                }
            }
        }

        std::cout << "palette size: " << palette.size() << std::endl;

        for (const auto& [iRow, rConflicts] : conflicts) {
            std::cout << iRow << " is in conflict with ";
            for (auto iColumn : rConflicts) std::cout << iColumn << " ";
            std::cout << std::endl;
        }
    }

    return 0;
}
