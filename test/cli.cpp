// --- External Includes ---
#include "mcgs/mcgs.hpp"
#include "parseMatrixMarket.hpp"

// --- STL Includes ---
#include <iostream> // std::cerr
#include <stdexcept>
#include <unordered_map> // std::unordered_map
#include <string> // std::string
#include <filesystem> // std::filesystem::path, std::filesystem::status, std::filesystem::file_type
#include <optional> // std::optional
#include <sstream> // std::stringstream
#include <cstring> // std::strlen
#include <fstream> // std::ifstream, std::ofstream


namespace mcgs {


const std::unordered_map<std::string,std::string> defaultArguments {
    {"-m", ""},
    {"-v", ""},
    {"-M", ""},
    {"-V", ""}
};


struct Arguments
{
    std::filesystem::path matrixInputPath;
    std::filesystem::path vectorInputPath;
    std::optional<std::filesystem::path> matrixOutputPath;
    std::optional<std::filesystem::path> vectorOutputPath;
}; // struct Arguments


struct InvalidArgument : public std::invalid_argument
{
    using std::invalid_argument::invalid_argument;
}; // struct InvalidArgument


void validateInputPath(const std::filesystem::path& rPath,
                       const std::string& rOption)
{
    const auto status = std::filesystem::status(rPath);
    switch (status.type()) {
        case std::filesystem::file_type::regular: break; // <== ok
        case std::filesystem::file_type::none: {
            std::stringstream message;
            message << "mcgscli: error: input file for option '"
                    << rOption << "' does not exist: " << rPath.string()
                    << '\n';
            throw InvalidArgument(message.str());
        }
        case std::filesystem::file_type::directory: {
            std::stringstream message;
            message << "mcgscli: error: input path for option '"
                    << rOption << "' is a directory.\n";
            throw InvalidArgument(message.str());
        }
        default:
            std::stringstream message;
            message << "mcgscli: error: unhandled file type for option '"
                    << rOption << "' at "
                    << rPath << "\n";
            throw InvalidArgument(message.str());
    }

    if ((status.permissions() & std::filesystem::perms::owner_read) == std::filesystem::perms::none) {
        std::stringstream message;
        message << "Error: missing read access to input file "
                << rPath
                << '\n';
        throw InvalidArgument(message.str());
    } // if !readPermission
}


void validateOutputPath(const std::filesystem::path& rPath,
                        const std::string& rOption,
                        bool allowOverwrite)
{
    const auto status = std::filesystem::status(rPath);
    switch (status.type()) {
        case std::filesystem::file_type::not_found: break; // <== ok
        case std::filesystem::file_type::regular: {
            if (!allowOverwrite) {
                std::stringstream message;
                message << "mcgscli: error: output path for option '"
                        << rOption << "' not exists: " << rPath.string()
                        << '\n';
                throw InvalidArgument(message.str());
            } else {
                break;
            }
        }
        case std::filesystem::file_type::directory: {
            std::stringstream message;
            message << "mcgscli: error: output path for option '"
                    << rOption << "' is a directory.\n";
            throw InvalidArgument(message.str());
        }
        default:
            std::stringstream message;
            message << "mcgscli: error: unhandled file type for option '"
                    << rOption << "' at "
                    << rPath << '\n';
            throw InvalidArgument(message.str());
    }

    if ((status.permissions() & std::filesystem::perms::owner_write) == std::filesystem::perms::none) {
        std::stringstream message;
        message << "Error: missing write access to output file "
                << rPath
                << '\n';
        throw InvalidArgument(message.str());
    } // if !writePermission
}


Arguments parseArguments(int argc, char const* const* argv)
{
    Arguments arguments;
    auto argMap = defaultArguments;

    // Parse optional argMap
    auto itArgument = argMap.end();
    for (int iArg=1; iArg<argc; ++iArg) {
        std::string arg = argv[iArg];
        if (!arg.empty() && arg.front() == '-') {
            // Parse a key
            if (itArgument == argMap.end()) {
                if ((itArgument = argMap.find(arg)) == argMap.end()) {
                    // The provided key does not exist in the argument map
                    throw InvalidArgument("Error: unrecognized option: " + arg + "\n");
                }
            } else {
                // No value was provided for the last key
                throw InvalidArgument("Error: missing argument for option " + itArgument->first + "\n");
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
                throw InvalidArgument("Error: missing option for argument " + arg + "\n");
            } // else (itArgument != argMap.end())
        } // else (!arg.empty() && arg.front() == '-')
    } // while (iArg < argc)

    // If the arg iterator was not reset, a value
    // was not provided for the last key.
    if (itArgument != argMap.end()) {
        throw InvalidArgument("Error: missing argument for option " + itArgument->first + "\n");
    }

    // Parse required arguments
    using PathString = std::filesystem::path::string_type;

    arguments.matrixInputPath = PathString(argMap["-m"]);
    validateInputPath(arguments.matrixInputPath, "-m");

    arguments.vectorInputPath = PathString(argMap["-v"]);
    validateInputPath(arguments.vectorInputPath, "-v");

    // Parse optional output arguments
    if (!argMap["-M"].empty()) {
        arguments.matrixOutputPath.emplace(PathString(argMap["-M"]));
        validateOutputPath(arguments.matrixOutputPath.value(), "-M", true);
    }

    if (!argMap["-V"].empty()) {
        arguments.vectorOutputPath.emplace(PathString(argMap["-V"]));
        validateOutputPath(arguments.vectorOutputPath.value(), "-V", true);
    }

    return arguments;
}


void print(const mcgs::CSRAdaptor<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value>& rMatrix,
           std::ostream* pStream = &std::cout)
{
    std::ostream& rStream = *pStream;
    rStream << std::fixed << std::setprecision(12) << std::scientific;
    rStream << "%%MatrixMarket matrix coordinate real general\n";
    rStream << rMatrix.rowCount << " " << rMatrix.columnCount << " " << rMatrix.entryCount << "\n";

    for (std::size_t iRow=0; iRow<rMatrix.rowCount; ++iRow) {
        const std::size_t iRowBegin = rMatrix.pRowExtents[iRow];
        const std::size_t iRowEnd = rMatrix.pRowExtents[iRow + 1];
        for (std::size_t iEntry=iRowBegin; iEntry<iRowEnd; ++iEntry) {
            const auto iColumn = rMatrix.pColumnIndices[iEntry];
            const auto entry = rMatrix.pEntries[iEntry];
            rStream << iRow + 1 << " " << iColumn + 1 << " " << entry << "\n";
        }
    }
}


int main(int argCount, char const* const* argValues)
{
    // Parse arguments
    const auto arguments = parseArguments(argCount, argValues);

    // Read input matrix and vector
    mcgs::TestCSRMatrix matrix;
    mcgs::TestDenseVector vector;

    {
        std::ifstream file(arguments.matrixInputPath);
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&matrix](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestCSRMatrix>) matrix = std::move(rValue);
        }, tmp);
    }

    {
        std::ifstream file(arguments.vectorInputPath);
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&vector](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestDenseVector>) vector = std::move(rValue);
        }, tmp);
    }

    if (matrix.columnCount != vector.size()) {
        std::stringstream message;
        message << "mcgscli: error: matrix(" << matrix.rowCount << "x" << matrix.columnCount << ") "
                << "- vector(" << vector.size() << "x1) size mismatch";
        throw InvalidArgument(message.str());
    }

    // Construct adaptor
    mcgs::CSRAdaptor<
        mcgs::TestCSRMatrix::Index,
        mcgs::TestCSRMatrix::Value
    > adaptor;
    adaptor.rowCount        = matrix.rowCount;
    adaptor.columnCount     = matrix.columnCount;
    adaptor.entryCount      = matrix.entryCount;
    adaptor.pRowExtents     = matrix.rowExtents.data();
    adaptor.pColumnIndices  = matrix.columnIndices.data();
    adaptor.pEntries        = matrix.entries.data();

    // Color
    std::vector<unsigned> colors(matrix.rowCount, std::numeric_limits<unsigned>::max());
    mcgs::ColorSettings<mcgs::TestCSRMatrix::Value> colorSettings;
    colorSettings.shrinkingFactor   = 256;
    colorSettings.maxStallCount     = 1e3;
    colorSettings.tolerance         = std::numeric_limits<mcgs::TestCSRMatrix::Value>::min();
    colorSettings.verbosity         = 1;

    mcgs::color(colors.data(), adaptor, colorSettings);

    // Partition and reorder
    auto pPartition = mcgs::makePartition<mcgs::TestCSRMatrix::Index>(colors.data(), adaptor.rowCount);
    auto pDummy = mcgs::reorder(matrix.rowCount, matrix.columnCount, matrix.entryCount,
                                matrix.rowExtents.data(), matrix.columnIndices.data(), matrix.entries.data(),
                                static_cast<double*>(nullptr), vector.data(),
                                pPartition);

    // Optional output
    if (arguments.matrixOutputPath.has_value()) {
        std::ofstream file(arguments.matrixOutputPath.value());
        print(adaptor, &file);
    }

    // Cleanup
    mcgs::destroyPartition(pPartition);
    mcgs::destroyPartition(pDummy);

    return MCGS_SUCCESS;
}


} // namespace mcgs


int main(int argCount, char const* const* argValues)
{
    try {
        return mcgs::main(argCount, argValues);
    } catch (mcgs::InvalidArgument& rException) {
        std::cerr << rException.what();
        return MCGS_FAILURE;
    } catch (std::exception& rException) {
        std::cerr << "mcgscli: mcgs terminated with the following exception message:\n"
                  << rException.what();
        return MCGS_FAILURE;
    } catch (...) {
        std::cerr << "mcgscli: mcgs terminated with an unknown exception.\n";
        return MCGS_FAILURE;
    }

    return MCGS_SUCCESS;
} // int main
