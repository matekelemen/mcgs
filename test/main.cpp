// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::Color, mcgs::Solve, mcgs::CSRAdaptor
#include "parseMatrixMarket.hpp" // mcgs::parseMatrixMarket

// --- STL Includes ---
#include <filesystem> // std::filesystem::path
#include <iostream> // std::cout, std::cerr
#include <fstream> // std::ifstream
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_set
#include <limits> // std::numeric_limits::max
#include <array> // std::array


void print(const mcgs::CSRAdaptor<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value>& rMatrix)
{
    std::ofstream file("matrix.mm");
    file << std::fixed << std::setprecision(12) << std::scientific;
    file << "%%MatrixMarket matrix coordinate real general\n";
    file << rMatrix.rowCount << " " << rMatrix.columnCount << " " << rMatrix.nonzeroCount << "\n";
    for (std::size_t iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        const std::size_t iRowBegin = rMatrix.pRowExtents[iRow];
        const std::size_t iRowEnd = rMatrix.pRowExtents[iRow + 1];
        for (std::size_t iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const auto iColumn = rMatrix.pColumnIndices[iEntry];
            const auto nonzero = rMatrix.pNonzeros[iEntry];
            file << iRow + 1 << " " << iColumn + 1 << " " << nonzero << "\n";
        }
    }
}


int main(int argc, const char* const * argv)
{
    if (argc != 3) {
        std::cerr << "Expecting exactly 2 arguments, but got " << argc << "\n";
        return MCGS_FAILURE;
    }

    // Read the input matrix and vector
    using Input = std::variant<mcgs::TestCSRMatrix,mcgs::TestDenseVector>;
    std::array<Input,2> input;

    for (std::size_t iArg=0ul; iArg<2ul; ++iArg) {
        const std::filesystem::path path = argv[iArg + 1];
        std::ifstream file(path);
        input[iArg] = mcgs::parseMatrixMarket(file);
    }

    mcgs::TestCSRMatrix* pMatrix = nullptr;
    mcgs::TestDenseVector* pVector = nullptr;
    for (auto& rInput : input) {
        std::visit([&pMatrix, &pVector](auto& rVariant){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rVariant)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestCSRMatrix>) {
                pMatrix = &rVariant;
            } else if constexpr (std::is_same_v<Type,mcgs::TestDenseVector>) {
                pVector = &rVariant;
            } else {
                std::cerr << "invalid input\n";
                return MCGS_FAILURE;
            }
        },
        rInput);
    }

    if (!pMatrix) {
        std::cerr << "missing input matrix\n";
        return MCGS_FAILURE;
    }

    if (!pVector) {
        std::cerr << "missing input vector\n";
        return MCGS_FAILURE;
    }

    // Construct a sparse CSR adaptor
    mcgs::CSRAdaptor<
        mcgs::TestCSRMatrix::Index,
        mcgs::TestCSRMatrix::Value
    > adaptor;
    adaptor.rowCount        = pMatrix->rowCount;
    adaptor.columnCount     = pMatrix->columnCount;
    adaptor.nonzeroCount    = pMatrix->nonzeroCount;
    adaptor.pRowExtents     = pMatrix->rowExtents.data();
    adaptor.pColumnIndices  = pMatrix->columnIndices.data();
    adaptor.pNonzeros       = pMatrix->nonzeros.data();

    // Allocate, initialize, define settings and color
    std::vector<unsigned> colors(pMatrix->columnCount, std::numeric_limits<unsigned>::max());

    {
        mcgs::ColorSettings settings;
        settings.verbosity = 1;
        settings.shrinkingFactor = 256;
        settings.maxStallCount = 1e4;
        mcgs::color(colors.data(), adaptor, settings);
    }

    // Check the coloring's correctness
    {
        std::unordered_map<
            mcgs::TestCSRMatrix::Index,
            std::vector<mcgs::TestCSRMatrix::Index>
        > conflicts;
        std::unordered_set<unsigned> palette;

        for (std::size_t iRow=0; iRow<pMatrix->columnCount; ++iRow) {
            const auto currentColor = colors[iRow];
            palette.insert(currentColor);

            const auto iRowBegin = pMatrix->rowExtents[iRow];
            const auto iRowEnd = pMatrix->rowExtents[iRow + 1];

            for (std::size_t iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
                const mcgs::TestCSRMatrix::Index iColumn = pMatrix->columnIndices[iEntry];
                const auto neighborColor = colors[iColumn];
                if (pMatrix->nonzeros[iEntry] && iRow != iColumn) {
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

    // Relax
    std::vector<mcgs::TestCSRMatrix::Value> solution(adaptor.columnCount, 0.0);
    {
        auto* pPartition = mcgs::makePartition(colors.data(), adaptor.columnCount);
        if (!pPartition) {
            std::cerr << "partitioning failed\n";
            return MCGS_FAILURE;
        }
        //mcgs::reorder(pMatrix->rowCount, pMatrix->columnCount, pMatrix->nonzeroCount,
        //              pMatrix->rowExtents.data(), pMatrix->columnIndices.data(), pMatrix->nonzeros.data(),
        //              pVector->data(),
        //              pPartition);
        //print(adaptor);

        mcgs::SolveSettings<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value> settings;
        settings.maxIterations = 1e1;
        settings.verbosity = 3;
        mcgs::solve(solution.data(), adaptor, pVector->data(), pPartition, settings);

        mcgs::destroyPartition<mcgs::TestCSRMatrix::Index>(pPartition);
    }

    return MCGS_SUCCESS;
}
