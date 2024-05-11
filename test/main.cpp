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
    if (argc != 2) {
        std::cerr << "Expecting exactly 1 argument, got " << argc << "\n";
        return MCGS_FAILURE;
    }

    // Read the input matrix
    mcgs::TestCSRMatrix matrix;
    {
        const std::filesystem::path matrixPath = argv[1];
        std::ifstream matrixFile(matrixPath);
        matrix = mcgs::parseMatrixMarket(matrixFile);
    }

    // Construct a sparse CSR adaptor
    mcgs::CSRAdaptor<
        mcgs::TestCSRMatrix::Index,
        mcgs::TestCSRMatrix::Value
    > adaptor;
    adaptor.rowCount = matrix.rowCount;
    adaptor.columnCount = matrix.columnCount;
    adaptor.nonzeroCount = matrix.nonzeroCount;
    adaptor.pRowExtents = matrix.rowExtents.data();
    adaptor.pColumnIndices = matrix.columnIndices.data();
    adaptor.pNonzeros = matrix.nonzeros.data();

    // Allocate, initialize, define settings and color
    std::vector<unsigned> colors(matrix.columnCount, std::numeric_limits<unsigned>::max());
    mcgs::ColorSettings settings;
    settings.verbosity = 3;
    //settings.shrinkingFactor = 5;
    mcgs::Color(adaptor, colors.data(), settings);

    // Check the coloring's correctness
    {
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
        if (!conflicts.empty()) {
            for (const auto& [iRow, rConflicts] : conflicts) {
                std::cerr << iRow << " is in conflict with ";
                for (auto iColumn : rConflicts) std::cerr << iColumn << " ";
                std::cerr << std::endl;
            }
            return MCGS_FAILURE;
        } // if conflicts
    } // destroy conflicts, palette

    return MCGS_SUCCESS;
}
