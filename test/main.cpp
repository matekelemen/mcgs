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
#include <chrono> // std::chrono::steady_clock, std::chrono::duration_cast


const std::unordered_map<std::string,std::string> defaultArguments {
    {"-s", "16"},
    {"-i", "1e2"}
};


struct Arguments
{
    std::filesystem::path matrixPath;
    std::filesystem::path vectorPath;
    int shrinkingFactor;
    std::size_t maxIterations;
}; // struct Arguments


Arguments parseArguments(int argc, char const* const* argv)
{
    Arguments arguments;
    auto argMap = defaultArguments;

    // Parse optional argMap
    auto itArgument = argMap.end();
    for (int iArg=3; iArg<argc; ++iArg) {
        std::string arg = argv[iArg];
        if (!arg.empty() && arg.front() == '-') {
            // Parse a key
            if (itArgument == argMap.end()) {
                if ((itArgument = argMap.find(arg)) == argMap.end()) {
                    // The provided key does not exist in the argument map
                    throw std::invalid_argument("Error: unrecognized option: " + arg + "\n");
                }
            } else {
                // No value was provided for the last key
                throw std::invalid_argument("Error: missing argument for option " + itArgument->first + "\n");
            } // else (itArgument == argMap.end())
        } else {
            // Parse a value
            if (itArgument != argMap.end()) {
                itArgument->second = arg;

                // Reset the arg iterator to indicate that a value was provided
                // for the current key.
                itArgument = argMap.end();
            } else {
                // No key was provided for this value
                throw std::invalid_argument("Error: missing option for argument " + arg + "\n");
            } // else (itArgument != argMap.end())
        } // else (!arg.empty() && arg.front() == '-')
    } // while (iArg < argc)

    // If the arg iterator was not reset, a value
    // was not provided for the last key.
    if (itArgument != argMap.end()) {
        throw std::invalid_argument("Error: missing argument for option " + itArgument->first + "\n");
    }

    // Parse required arguments
    if (argc < 2) {
        throw std::invalid_argument("Error: missing input matrix\n");
    } else if (argc < 3) {
        throw std::invalid_argument("Error: missing input vector\n");
    }
    using PathString = std::filesystem::path::string_type;
    arguments.matrixPath = PathString(argv[1]);
    arguments.vectorPath = PathString(argv[2]);

    // Validate input path
    for (const auto& rPath : {arguments.matrixPath, arguments.vectorPath}) {
        const auto status = std::filesystem::status(rPath);
        switch (status.type()) {
            case std::filesystem::file_type::regular: break; // <== ok
            case std::filesystem::file_type::none: throw std::invalid_argument(
                "Error: input file does not exist: " + rPath.string() + "\n");
            case std::filesystem::file_type::directory: throw std::invalid_argument(
                "Error: provided input path is a directory: " + rPath.string() + "\n");
            default: throw std::invalid_argument(
                "Error: input is not a file: " + rPath.string() + "\n");
            } // switch status

        if ((status.permissions() & std::filesystem::perms::owner_read) == std::filesystem::perms::none) {
            throw std::invalid_argument(
                "Error: missing read access to input file " + rPath.string() + "\n");
        } // if !readPermission
    }

    // Convert and validate shrinking factor
    char* itEnd = nullptr;
    const std::string& rShrinkingFactorString = argMap["-s"];
    const int shrinkingFactor = std::strtol(rShrinkingFactorString.data(), &itEnd, 0);
    if (itEnd < rShrinkingFactorString.data() ||
        static_cast<std::size_t>(std::distance(rShrinkingFactorString.data(), static_cast<const char*>(itEnd))) != rShrinkingFactorString.size()) {
        throw std::invalid_argument(
            "Error: invalid shrinking factor: " + rShrinkingFactorString + "\n");
    } else {
        arguments.shrinkingFactor = static_cast<std::size_t>(shrinkingFactor);
    }

    // Convert and validate max iterations
    itEnd = nullptr;
    const std::string& rMaxIterationsString = argMap["-i"];
    const std::size_t maxIterations = std::strtod(rMaxIterationsString.data(), &itEnd);
    if (itEnd < rMaxIterationsString.data() ||
        static_cast<std::size_t>(std::distance(rMaxIterationsString.data(), static_cast<const char*>(itEnd))) != rMaxIterationsString.size()) {
        throw std::invalid_argument(
            "Error: invalid max iterations: " + rMaxIterationsString + "\n");
    } else {
        arguments.maxIterations = static_cast<std::size_t>(maxIterations);
    }

    return arguments;
}


class ScopedTimer
{
public:
    ScopedTimer(std::string&& rTag)
        : _begin(std::chrono::steady_clock::now()),
          _tag(std::move(rTag))
    {}

    ~ScopedTimer()
    {
        const auto end = std::chrono::steady_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - _begin);
        std::cout << _tag << " took " << elapsed.count() << "[ms]\n";
    }

private:
    std::chrono::steady_clock::time_point _begin;
    std::string _tag;
};


#define MCGS_SCOPED_TIMER(TAG) [[maybe_unused]] ScopedTimer MCGS_SCOPED_TIMER_INSTANCE(TAG)


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
    // --- Reading ---
    const auto arguments = parseArguments(argc, argv);

    // Read the input matrix and vector
    mcgs::TestCSRMatrix matrix;
    mcgs::TestDenseVector vector;

    {
        MCGS_SCOPED_TIMER("Reading " + arguments.matrixPath.string());
        std::ifstream file(arguments.matrixPath);
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&matrix](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestCSRMatrix>) matrix = std::move(rValue);
        }, tmp);
    }

    {
        MCGS_SCOPED_TIMER("Reading " + arguments.vectorPath.string());
        std::ifstream file(arguments.vectorPath);
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&vector](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestDenseVector>) vector = std::move(rValue);
        }, tmp);
    }

    // Construct a sparse CSR adaptor
    mcgs::CSRAdaptor<
        mcgs::TestCSRMatrix::Index,
        mcgs::TestCSRMatrix::Value
    > adaptor;
    adaptor.rowCount        = matrix.rowCount;
    adaptor.columnCount     = matrix.columnCount;
    adaptor.nonzeroCount    = matrix.nonzeroCount;
    adaptor.pRowExtents     = matrix.rowExtents.data();
    adaptor.pColumnIndices  = matrix.columnIndices.data();
    adaptor.pNonzeros       = matrix.nonzeros.data();

    // --- Coloring ---
    std::vector<unsigned> colors(matrix.columnCount, std::numeric_limits<unsigned>::max());
    mcgs::ColorSettings<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value> colorSettings;

    {
        MCGS_SCOPED_TIMER("coloring");
        colorSettings.verbosity = 2;
        colorSettings.shrinkingFactor = arguments.shrinkingFactor;
        colorSettings.maxStallCount = 1e4;
        colorSettings.tolerance = 1e-10;
        mcgs::color(colors.data(), adaptor, colorSettings);
    }

    // Check the coloring's correctness
    {
        MCGS_SCOPED_TIMER("coloring validation");
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
                if (colorSettings.tolerance < std::abs(matrix.nonzeros[iEntry]) && iRow != iColumn) {
                    if (neighborColor == currentColor) {
                        conflicts.emplace(iRow, std::vector<mcgs::TestCSRMatrix::Index> {})
                            .first->second.push_back(iColumn);
                    }
                }
            }
        }

        std::cout << "palette size: " << palette.size() << std::endl;
        if (!conflicts.empty()) {
            //for (const auto& [iRow, rConflicts] : conflicts) {
            //    std::cerr << iRow << " is in conflict with ";
            //    for (auto iColumn : rConflicts) std::cerr << iColumn << " ";
            //    std::cerr << std::endl;
            //}
            std::cerr << "coloring failed; " << conflicts.size() << " rows are in conflict\n";
            return MCGS_FAILURE;
        } // if conflicts
    } // destroy conflicts, palette

    // --- Partitioning ---
    mcgs::Partition<mcgs::TestCSRMatrix::Index>* pPartition = nullptr;

    {
        MCGS_SCOPED_TIMER("partitioning");
        pPartition = mcgs::makePartition(colors.data(), adaptor.columnCount);
        if (!pPartition) {
            std::cerr << "partitioning failed\n";
            return MCGS_FAILURE;
        }
    }

    // --- Relaxation ---
    std::vector<mcgs::TestCSRMatrix::Value> solution(adaptor.columnCount);
    mcgs::SolveSettings<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value> settings;
    settings.maxIterations = arguments.maxIterations;
    settings.verbosity = 1;

    std::vector<mcgs::TestCSRMatrix::Value> buffer(adaptor.columnCount);
    const auto initialResidual = mcgs::residual(adaptor, solution.data(), vector.data());

    // Serial relaxation
    std::fill(solution.begin(), solution.end(), 0.0);
    {
        MCGS_SCOPED_TIMER("serial relaxation");
        settings.parallelization = mcgs::Parallelization::None;
        if (mcgs::solve(solution.data(), adaptor, vector.data(), settings) != MCGS_SUCCESS) {
            std::cerr << "serial relaxation failed\n";
            return MCGS_FAILURE;
        }
    }
    std::cout << "residual " << mcgs::residual(adaptor, solution.data(), vector.data()) / initialResidual << "\n";

    // Parallel relaxation
    mcgs::Partition<mcgs::TestCSRMatrix::Index>* pReorderedPartition = nullptr;
    {
        MCGS_SCOPED_TIMER("reordering");
        pReorderedPartition = mcgs::reorder(matrix.rowCount, matrix.columnCount, matrix.nonzeroCount,
                                            matrix.rowExtents.data(), matrix.columnIndices.data(), matrix.nonzeros.data(),
                                            vector.data(),
                                            pPartition);
        if (!pReorderedPartition) {
            std::cerr << "reordering failed\n";
            return MCGS_FAILURE;
        }
    }

    std::fill(solution.begin(), solution.end(), 0.0);
    {
        MCGS_SCOPED_TIMER("row-wise parallel relaxation");
        settings.parallelization = mcgs::Parallelization::RowWise;
        if (mcgs::solve(solution.data(), adaptor, vector.data(), pReorderedPartition, settings) != MCGS_SUCCESS) {
            std::cerr << "parallel relaxation failed\n";
            return MCGS_FAILURE;
        }
    }
    std::cout << "residual " << mcgs::residual(adaptor, solution.data(), vector.data()) / initialResidual << "\n";

    std::fill(solution.begin(), solution.end(), 0.0);
    {
        MCGS_SCOPED_TIMER("entrywise parallel relaxation");
        settings.parallelization = mcgs::Parallelization::EntryWise;
        if (mcgs::solve(solution.data(), adaptor, vector.data(), pReorderedPartition, settings) != MCGS_SUCCESS) {
            std::cerr << "parallel relaxation failed\n";
            return MCGS_FAILURE;
        }
    }
    std::cout << "residual " << mcgs::residual(adaptor, solution.data(), vector.data()) / initialResidual << "\n";

    {
        MCGS_SCOPED_TIMER("undo reordering");
        mcgs::revertReorder(solution.data(), pPartition);
    }

    // Cleanup
    mcgs::destroyPartition<mcgs::TestCSRMatrix::Index>(pReorderedPartition);
    mcgs::destroyPartition<mcgs::TestCSRMatrix::Index>(pPartition);

    return MCGS_SUCCESS;
}
