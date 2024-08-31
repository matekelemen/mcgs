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
#include <sstream> // std::basic_stringstream
#include <cstring> // std::strlen
#include <fstream> // std::ifstream, std::ofstream


namespace mcgs {


using PathString = std::filesystem::path::string_type;
using PathStream = std::basic_stringstream<PathString::value_type>;
using ArgumentMap = std::unordered_map<
    std::string,    ///< option name
    PathString      ///< option value
>;


PathString makePathString(const char* pPath)
{
    PathString output(pPath, pPath + std::strlen(pPath));
    return output;
}


const ArgumentMap defaultArguments {
    {"-a", makePathString("")},
    {"-x", makePathString("")},
    {"-b", makePathString("")},
    {"-A", makePathString("")},
    {"-X", makePathString("")},
    {"-s", makePathString("16")},
    {"-i", makePathString("1e2")},
    {"-r", makePathString("1")}
};


void help(std::ostream& rStream)
{
}


struct Arguments
{
    std::optional<std::filesystem::path> lhsInputPath;
    std::optional<std::filesystem::path> solutionInputPath;
    std::optional<std::filesystem::path> rhsInputPath;
    std::optional<std::filesystem::path> lhsOutputPath;
    std::optional<std::filesystem::path> solutionOutputPath;
    int shrinkingFactor;
    int maxIterations;
    double relaxation;
}; // struct Arguments


struct InvalidArgument : public std::invalid_argument
{
    using std::invalid_argument::invalid_argument;
}; // struct InvalidArgument


struct RuntimeError : public std::runtime_error
{
    using std::runtime_error::runtime_error;
}; // struct RuntimeError


InvalidArgument makeInvalidArgument(const PathString& rMessage)
{
    return InvalidArgument(std::string(rMessage.begin(), rMessage.end()));
}


void parseInputPath(std::optional<std::filesystem::path>& rPath,
                    const ArgumentMap& rArguments,
                    const ArgumentMap::key_type& rOption)
{
    std::filesystem::path path = rArguments.at(rOption);

    if (!path.empty()) {
        const auto status = std::filesystem::status(path);
        switch (status.type()) {
            case std::filesystem::file_type::regular: break; // <== ok
            case std::filesystem::file_type::none: {
                PathStream message;
                message << "mcgscli: error: input file for option '"
                        << PathString(rOption.begin(), rOption.end()) << "' does not exist: " << path
                        << '\n';
                throw makeInvalidArgument(message.str());
            }
            case std::filesystem::file_type::directory: {
                PathStream message;
                message << "mcgscli: error: input path for option '"
                        << PathString(rOption.begin(), rOption.end()) << "' is a directory.\n";
                throw makeInvalidArgument(message.str());
            }
            default:
                PathStream message;
                message << "mcgscli: error: unhandled file type for option '"
                        << PathString(rOption.begin(), rOption.end()) << "' at "
                        << path << "\n";
                throw makeInvalidArgument(message.str());
        }

        if ((status.permissions() & std::filesystem::perms::owner_read) == std::filesystem::perms::none) {
            PathStream message;
            message << "Error: missing read access to input file "
                    << path
                    << '\n';
            throw makeInvalidArgument(message.str());
        } // if !readPermission

        rPath.emplace(std::move(path));
    } // if path
}


void parseOutputPath(std::optional<std::filesystem::path>& rPath,
                     const ArgumentMap& rArguments,
                     const ArgumentMap::key_type& rOption,
                     bool allowOverwrite)
{
    std::filesystem::path path(rArguments.at(rOption).begin(), rArguments.at(rOption).end());

    if (!path.empty()) {
        const auto status = std::filesystem::status(path);
        switch (status.type()) {
            case std::filesystem::file_type::not_found: break; // <== ok
            case std::filesystem::file_type::regular: {
                if (!allowOverwrite) {
                    PathStream message;
                    message << "mcgscli: error: output path for option '"
                            << PathString(rOption.begin(), rOption.end()) << "' not exists: " << path
                            << '\n';
                    throw makeInvalidArgument(message.str());
                } else {
                    break;
                }
            }
            case std::filesystem::file_type::directory: {
                PathStream message;
                message << "mcgscli: error: output path for option '"
                        << PathString(rOption.begin(), rOption.end()) << "' is a directory.\n";
                throw makeInvalidArgument(message.str());
            }
            default:
                PathStream message;
                message << "mcgscli: error: unhandled file type for option '"
                        << PathString(rOption.begin(), rOption.end()) << "' at "
                        << path << '\n';
                throw makeInvalidArgument(message.str());
        }

        if ((status.permissions() & std::filesystem::perms::owner_write) == std::filesystem::perms::none) {
            PathStream message;
            message << "Error: missing write access to output file "
                    << path
                    << '\n';
            throw makeInvalidArgument(message.str());
        } // if !writePermission

        rPath.emplace(std::move(path));
    } // if path
}


template <class T>
void parseNumericArgument(T& rValue,
                          const ArgumentMap& rArguments,
                          const ArgumentMap::key_type& rOption)
{
    PathStream stream;
    stream << rArguments.at(rOption);
    try {
        long double tmp;
        stream >> tmp;
        rValue = tmp;
    } catch (std::exception& rException) {
        stream.clear();
        stream << "mcgscli: error: invalid argument for option '"
               << PathString(rOption.begin(), rOption.end()) << "': "
               << rArguments.at(rOption)
               << "\n";
        throw makeInvalidArgument(stream.str());
    }
}


Arguments parseArguments(int argc, char const* const* argv)
{
    Arguments arguments;
    auto argMap = defaultArguments;

    // Parse optional argMap
    auto itArgument = argMap.end();
    for (int iArg=1; iArg<argc; ++iArg) {
        PathString arg(argv[iArg], argv[iArg] + std::strlen(argv[iArg]));
        if (!arg.empty() && arg.front() == '-') {
            // Parse a key
            ArgumentMap::key_type key(arg.begin(), arg.end());
            if (itArgument == argMap.end()) {
                if ((itArgument = argMap.find(key)) == argMap.end()) {
                    // Special case if '-h' or '--help'
                    if (key == "-h" || key == "--help") {
                        help(std::cout);
                    } else {
                        help(std::cerr);

                        // The provided key does not exist in the argument map
                        PathStream stream;
                        stream << "mcgscli: error: unrecognized option: "
                            << arg
                            << '\n';
                        throw makeInvalidArgument(stream.str());
                    }
                }
            } else {
                // No value was provided for the last key
                PathStream stream;
                stream << "mcgscli: error: missing argument for option "
                       << PathString(itArgument->first.begin(), itArgument->first.end())
                       << '\n';
                throw makeInvalidArgument(stream.str());
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
                PathStream stream;
                stream << "mcgscli: error: missing option for argument "
                       << arg
                       << '\n';
                throw makeInvalidArgument(stream.str());
            } // else (itArgument != argMap.end())
        } // else (!arg.empty() && arg.front() == '-')
    } // while (iArg < argc)

    // If the arg iterator was not reset, a value
    // was not provided for the last key.
    if (itArgument != argMap.end()) {
        PathStream stream;
        stream << "mcgscli: error: missing argument for option "
               << PathString(itArgument->first.begin(), itArgument->first.end())
               << '\n';
        throw makeInvalidArgument(stream.str());
    }

    // Parse input paths.
    parseInputPath(arguments.lhsInputPath, argMap, "-a");
    parseInputPath(arguments.solutionInputPath, argMap, "-x");
    parseInputPath(arguments.rhsInputPath, argMap, "-b");
    
    // Parse output paths.
    parseOutputPath(arguments.lhsOutputPath, argMap, "-A", true);
    parseOutputPath(arguments.solutionOutputPath, argMap, "-X", true);

    // Parse other arguments.
    parseNumericArgument(arguments.shrinkingFactor, argMap, "-s");
    parseNumericArgument(arguments.maxIterations, argMap, "-i");
    parseNumericArgument(arguments.relaxation, argMap, "-r");

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


void print(const mcgs::TestDenseVector& rVector,
           std::ostream* pStream = &std::cout)
{
    std::ostream& rStream = *pStream;
    rStream << std::fixed << std::setprecision(12) << std::scientific
            << "%%MatrixMarket matrix array real general\n"
            << rVector.size() << " 1\n";
    for (const auto item : rVector) rStream << item << '\n';
}


int main(int argCount, char const* const* argValues)
{
    // Parse arguments
    const auto arguments = parseArguments(argCount, argValues);

    // Read input matrix and vectors
    std::optional<mcgs::TestCSRMatrix> lhs;
    std::optional<mcgs::TestDenseVector> solution, rhs;

    if (arguments.lhsInputPath.has_value()) {
        std::ifstream file(arguments.lhsInputPath.value());
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&lhs](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestCSRMatrix>) lhs = std::move(rValue);
        }, tmp);
    }

    if (arguments.rhsInputPath.has_value()) {
        std::ifstream file(arguments.rhsInputPath.value());
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&rhs](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestDenseVector>) rhs = std::move(rValue);
        }, tmp);
    }

    if (arguments.solutionInputPath.has_value()) {
        std::ifstream file(arguments.solutionInputPath.value());
        auto tmp = mcgs::parseMatrixMarket(file);
        std::visit([&solution](auto& rValue){
            using Type = std::remove_const_t<std::remove_reference_t<decltype(rValue)>>;
            if constexpr (std::is_same_v<Type,mcgs::TestDenseVector>) solution = std::move(rValue);
        }, tmp);
    }

    // Sanity checks
    if (lhs.has_value() && rhs.has_value()) {
        if (lhs.value().columnCount != rhs.value().size()) {
            PathStream message;
            message << "mcgscli: error: matrix(" << lhs.value().rowCount << "x" << lhs.value().columnCount << ") "
                    << "- vector(" << rhs.value().size() << "x1) size mismatch";
            throw makeInvalidArgument(message.str());
        }
    }

    if (lhs.has_value() && solution.has_value()) {
        if (lhs.value().columnCount != solution.value().size()) {
            PathStream message;
            message << "mcgscli: error: matrix(" << lhs.value().rowCount << "x" << lhs.value().columnCount << ") "
                    << "- vector(" << solution.value().size() << "x1) size mismatch";
            throw makeInvalidArgument(message.str());
        }
    }

    if (solution.has_value() && rhs.has_value()) {
        if (solution.value().size() != rhs.value().size()) {
            PathStream message;
            message << "mcgscli: error: solution (" << solution.value().size() << "x1) "
                    << "- right hand side (" << rhs.value().size() << "x1) size mismatch";
            throw makeInvalidArgument(message.str());
        }
    }

    if (lhs.has_value()) {
        // Construct adaptor
        mcgs::CSRAdaptor<
            mcgs::TestCSRMatrix::Index,
            mcgs::TestCSRMatrix::Value
        > adaptor;
        adaptor.rowCount        = lhs.value().rowCount;
        adaptor.columnCount     = lhs.value().columnCount;
        adaptor.entryCount      = lhs.value().entryCount;
        adaptor.pRowExtents     = lhs.value().rowExtents.data();
        adaptor.pColumnIndices  = lhs.value().columnIndices.data();
        adaptor.pEntries        = lhs.value().entries.data();

        // Color
        std::vector<unsigned> colors(lhs.value().rowCount, std::numeric_limits<unsigned>::max());
        mcgs::ColorSettings<mcgs::TestCSRMatrix::Value> colorSettings;
        colorSettings.shrinkingFactor   = arguments.shrinkingFactor;
        colorSettings.maxStallCount     = 1e3;
        colorSettings.tolerance         = std::numeric_limits<mcgs::TestCSRMatrix::Value>::min();
        colorSettings.verbosity         = 1;

        mcgs::color(colors.data(), adaptor, colorSettings);

        // Partition and reorder
        double* pSolution = solution.has_value() ? solution.value().data() : nullptr;
        double* pRhs = rhs.has_value() ? rhs.value().data() : nullptr;

        auto pPartition = mcgs::makePartition<mcgs::TestCSRMatrix::Index>(colors.data(), adaptor.rowCount);
        auto pReorderedPartition = mcgs::reorder(lhs.value().rowCount, lhs.value().columnCount, lhs.value().entryCount,
                                                 lhs.value().rowExtents.data(), lhs.value().columnIndices.data(), lhs.value().entries.data(),
                                                 pSolution, pRhs,
                                                 pPartition);

        // Optional output
        if (arguments.lhsOutputPath.has_value()) {
            std::ofstream file(arguments.lhsOutputPath.value());
            print(adaptor, &file);
        }

        // Smoothing
        if (rhs.has_value()) {
            if (!solution.has_value()) {
                solution.emplace(mcgs::TestDenseVector(lhs.value().rowCount, 0.0));
            }

            mcgs::SolveSettings<mcgs::TestCSRMatrix::Index,mcgs::TestCSRMatrix::Value> settings;
            settings.maxIterations = arguments.maxIterations;
            settings.relaxation = arguments.relaxation;
            settings.parallelization = mcgs::Parallelization::EntryWise;
            settings.verbosity = 1;
            if (mcgs::solve(solution.value().data(),
                            adaptor,
                            rhs.value().data(),
                            pReorderedPartition,
                            settings) != MCGS_SUCCESS) {
                throw RuntimeError("mcgscli: error: smoothing failed\n");
            }

            if (arguments.solutionOutputPath.has_value()) {
                std::ofstream file(arguments.solutionOutputPath.value());
                print(solution.value(), &file);
            }
        }

        // Cleanup
        mcgs::destroyPartition(pPartition);
        mcgs::destroyPartition(pReorderedPartition);
    } /*if lhs.has_value()*/

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
    } catch(mcgs::RuntimeError& rException) {
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
