// --- Internal Includes ---
#include "parseMatrixMarket.hpp"

// --- STL Includes ---
#include <istream> // std::istream
#include <type_traits>
#include <variant> // std::variant
#include <regex> // std::regex
#include <stdexcept> // std::runtime_error
#include <string> // std::string, std::to_string
#include <limits> // std::numeric_limits::max

namespace mcgs {


std::variant<
    TestCSRMatrix,
    TestDenseVector
> parseMatrixMarketHeader(std::istream& rStream)
{
    std::variant<TestCSRMatrix,TestDenseVector> output;

    std::regex formatPattern(R"(^%%MatrixMarket (\w+) (\w+) (.*)?)");
    std::regex qualifierPattern(R"(\w+)");
    std::size_t iLine = 0ul;
    std::string buffer(0x400, '\0');

    while (rStream.peek() == '%' && !rStream.bad() && !rStream.fail()) /*comment line begins with a '%'*/ {
        // Read the input stream until the buffer is filled
        // or a newline character is encountered.
        rStream.getline(buffer.data(), buffer.size());

        // If the fail bit is set, the buffer got filled before
        // a newline character was found in the stream. This means
        // that the input file is invalid.
        if (rStream.fail()) throw std::runtime_error("Error: failed to read line " + std::to_string(iLine));

        // The first line must contain format properties.
        std::match_results<const std::string::value_type*> match;
        if (std::regex_match(buffer.data(), match, formatPattern)) {
            if (iLine == 0ul) {
                if (match.size() < 2) throw std::runtime_error("Error: too few matches in line " + std::to_string(iLine));

                // Parse object type (matrix, vector, etc.)
                const std::string objectName = match.str(1);
                if (objectName != "matrix") throw std::runtime_error("Error: invalid object type: " + objectName);

                // Parse format type (coordinate or array)
                const std::string formatName = match.str(2);
                if (formatName == "coordinate") {
                    output = TestCSRMatrix();
                } else if (formatName == "array") {
                    output = TestDenseVector();
                } else {
                    throw std::runtime_error("Error: invalid format type: " + formatName);
                }

                // Loop over optional qualifiers (value type, symmetry)
                if (3 < match.size()) {
                    const std::string qualifiers = match.str(3);
                    for (auto itQualifier = std::sregex_iterator(qualifiers.begin(),
                                                                 qualifiers.end(),
                                                                 qualifierPattern);
                            itQualifier != std::sregex_iterator();
                            ++itQualifier) {
                        const std::string qualifier = itQualifier->str();
                        if (qualifier == "real") {
                            // ok
                        } else if (qualifier == "integer") {
                            throw std::runtime_error("Error: integer matrices are unsupported");
                        } else if (qualifier == "complex") {
                            throw std::runtime_error("Error: complex matrices are unsupported");
                        } else if (qualifier == "pattern") {
                            throw std::runtime_error("Pattern?");
                        } else if (qualifier == "general") {
                            // ok
                        } else if (qualifier == "symmetric") {
                            throw std::runtime_error("Error: symmetric matrices are unsupported");
                        } else if (qualifier == "skew-symmetric") {
                            throw std::runtime_error("Error: skew-symmetric matrices are unsupported");
                        } else if (qualifier == "hermitian") {
                            throw std::runtime_error("Error: hermitian matrices are unsupported");
                        } else {
                            throw std::runtime_error("Error: invalid qualifier: " + qualifier);
                        }
                    }
                }
            } // if iLine == 0

        } else if (iLine == 0ul) {
            throw std::runtime_error("Error: the first line of the input must begin with '%%MatrixMarket' and define the matrix format");
        }

        ++iLine;
    }

    long long rows = 0ul, columns = 0ul, nonzeros = 0ul;

    // Parse number of rows
    rStream >> rows;
    if (rStream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of rows in the input matrix");
    } else if (rows < 0ll) {
        throw std::runtime_error("Error: negative number of rows in input");
    } else {
        std::visit([rows](auto& rVariant){
                using Value = std::remove_reference_t<decltype(rVariant)>;
                if constexpr (std::is_same_v<Value,TestCSRMatrix>) {
                    rVariant.rowCount = rows;
                } else if constexpr (std::is_same_v<Value,TestDenseVector>) {
                    rVariant.reserve(rows);
                }
            }, output);
    }

    // Parse number of columns
    rStream >> columns;
    if (rStream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of columns in the input matrix");
    } else if (columns < 0ll) {
        throw std::runtime_error("Error: negative number of columns in input");
    } else {
        std::visit([columns](auto& rVariant){
                using Value = std::remove_reference_t<decltype(rVariant)>;
                if constexpr (std::is_same_v<Value,TestCSRMatrix>) {
                    rVariant.columnCount = columns;
                } else if constexpr (std::is_same_v<Value,TestDenseVector>) {
                    if (columns != 1) throw std::runtime_error("Error: expecting a vector, but the column size is " + std::to_string(columns));
                }
            }, output);
    }

    // Parse number of nonzeros
    rStream >> nonzeros;
    if (rStream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of nonzeros in the input matrix\n");
    } else if (nonzeros < 0ll) {
        throw std::runtime_error("Error: negative number of nonzeros in input");
    } else {
        std::visit([nonzeros](auto& rVariant){
                using Value = std::remove_reference_t<decltype(rVariant)>;
                if constexpr (std::is_same_v<Value,TestCSRMatrix>) {
                    rVariant.nonzeroCount = nonzeros;
                }
            }, output);
    }

    // Ignore the rest of the line
    rStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::visit([](auto& rVariant){
            using Value = std::remove_reference_t<decltype(rVariant)>;
            if constexpr (std::is_same_v<Value,TestCSRMatrix>) {
                rVariant.rowExtents.reserve(rVariant.rowCount + 1);
                rVariant.columnIndices.reserve(rVariant.nonzeroCount);
                rVariant.nonzeros.reserve(rVariant.nonzeroCount);
            }
        }, output);

    return output;
}


bool parseSparseDataLine(std::istream& rStream, TestCSRMatrix& rMatrix)
{
    constexpr std::streamsize ignoreSize = std::numeric_limits<std::streamsize>::max();
    TestCSRMatrix::Index iRow = 0, iColumn = 0;
    TestCSRMatrix::Value value = 0.0;

    // Read row index
    rStream >> iRow;
    if (rStream.eof() || rStream.bad()) {
        rStream.clear();
        rStream.ignore(ignoreSize, '\n');
        return false;
    }

    // Read column index
    rStream >> iColumn;
    if (rStream.eof() || rStream.bad()) {
        rStream.clear();
        rStream.ignore(ignoreSize, '\n');
        return false;
    }

    // Read value if requested
    rStream >> value;
    if (rStream.eof() || rStream.bad()) {
        rStream.clear();
        rStream.ignore(ignoreSize, '\n');
        return false;
    }

    rStream.ignore(ignoreSize, '\n');

    // Convert (row,column) indices to 0-based indices
    --iRow;
    --iColumn;

    // Record entry in the output matrix
    for (auto i=rMatrix.rowExtents.size(); i<=iRow; ++i) rMatrix.rowExtents.push_back(rMatrix.nonzeros.size());
    rMatrix.columnIndices.push_back(iColumn);
    rMatrix.nonzeros.push_back(value);

    return true;
}


bool parseDenseDataLine(std::istream& rStream, TestDenseVector& rVector)
{
    constexpr std::streamsize ignoreSize = std::numeric_limits<std::streamsize>::max();

    rVector.emplace_back();
    rStream >> rVector.back();
    if (rStream.eof() || rStream.bad()) {
        rStream.clear();
        rStream.ignore(ignoreSize, '\n');
        rVector.pop_back();
        return false;
    }

    return true;
}


std::variant<
    TestCSRMatrix,
    TestDenseVector
> parseMatrixMarket(std::istream& rStream)
{
    auto input = parseMatrixMarketHeader(rStream);

    std::visit([&rStream](auto& rVariant){
        using Type = std::remove_reference_t<decltype(rVariant)>;
        if constexpr (std::is_same_v<Type,TestCSRMatrix>) {
            while (parseSparseDataLine(rStream, rVariant)) {}
            rVariant.rowExtents.push_back(rVariant.nonzeros.size());
        } else if constexpr (std::is_same_v<Type,TestDenseVector>) {
            while (parseDenseDataLine(rStream, rVariant)) {}
        }
    }, input);

    return input;
}


} // namespace mcgs
