// --- Internal Includes ---
#include "parseMatrixMarket.hpp"

// --- STL Includes ---
#include <istream> // std::istream
#include <regex> // std::regex
#include <stdexcept> // std::runtime_error
#include <string> // std::string, std::to_string
#include <limits> // std::numeric_limits::max

namespace mcgs {


TestCSRMatrix parseMatrixMarketHeader(std::istream& r_stream)
{
    TestCSRMatrix output;

    std::regex formatPattern(R"(^%%MatrixMarket (\w+) (\w+) (.*)?)");
    std::regex qualifierPattern(R"(\w+)");
    std::size_t i_line = 0ul;
    std::string buffer(0x400, '\0');

    while (r_stream.peek() == '%' && !r_stream.bad() && !r_stream.fail()) /*comment line begins with a '%'*/ {
        // Read the input stream until the buffer is filled
        // or a newline character is encountered.
        r_stream.getline(buffer.data(), buffer.size());

        // If the fail bit is set, the buffer got filled before
        // a newline character was found in the stream. This means
        // that the input file is invalid.
        if (r_stream.fail()) throw std::runtime_error("Error: failed to read line " + std::to_string(i_line));

        // The first line must contain format properties.
        std::match_results<const std::string::value_type*> match;
        if (std::regex_match(buffer.data(), match, formatPattern)) {
            if (i_line == 0ul) {
                if (match.size() < 2) throw std::runtime_error("Error: too few matches in line " + std::to_string(i_line));

                // Parse object type (matrix, vector, etc.)
                const std::string objectName = match.str(1);
                if (objectName != "matrix") throw std::runtime_error("Error: invlaid object type: " + objectName);

                // Parse format type (coordinate or array)
                const std::string formatName = match.str(2);
                if (formatName != "coordinate") throw std::runtime_error("Error: invalid format type: " + formatName);

                // Loop over optional qualifiers (value type, symmetry)
                if (3 < match.size()) {
                    const std::string qualifiers = match.str(3);
                    for (auto it_qualifier = std::sregex_iterator(qualifiers.begin(),
                                                                  qualifiers.end(),
                                                                  qualifierPattern);
                            it_qualifier != std::sregex_iterator();
                            ++it_qualifier) {
                        const std::string qualifier = it_qualifier->str();
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
            } // if i_line == 0

        } else if (i_line == 0ul) {
            throw std::runtime_error("Error: the first line of the input must begin with '%%MatrixMarket' and define the matrix format");
        }

        ++i_line;
    }

    long long rows = 0ul, columns = 0ul, nonzeros = 0ul;

    // Parse number of rows
    r_stream >> rows;
    if (r_stream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of rows in the input matrix");
    } else if (rows < 0ll) {
        throw std::runtime_error("Error: negative number of rows in input");
    } else {
        output.rowCount = static_cast<TestCSRMatrix::Index>(rows);
    }

    // Parse number of columns
    r_stream >> columns;
    if (r_stream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of columns in the input matrix");
    } else if (columns < 0ll) {
        throw std::runtime_error("Error: negative number of columns in input");
    } else {
        output.columnCount = static_cast<std::size_t>(columns);
    }

    // Parse number of nonzeros
    r_stream >> nonzeros;
    if (r_stream.fail()) {
        throw std::runtime_error("Error: failed to parse the number of nonzeros in the input matrix\n");
    } else if (nonzeros < 0ll) {
        throw std::runtime_error("Error: negative number of nonzeros in input");
    } else {
        output.nonzeroCount = static_cast<std::size_t>(nonzeros);
    }

    // Ignore the rest of the line
    r_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    output.rowExtents.reserve(output.rowCount + 1);
    output.rowExtents.push_back(0);
    output.columnIndices.reserve(output.nonzeroCount);
    output.nonzeros.reserve(output.nonzeroCount);
    return output;
}


bool parseSparseDataLine(std::istream& r_stream, TestCSRMatrix& r_matrix)
{
    constexpr std::streamsize ignoreSize = std::numeric_limits<std::streamsize>::max();
    TestCSRMatrix::Index i_row = 0, i_column = 0;
    TestCSRMatrix::Value value = 0.0;

    // Read row index
    r_stream >> i_row;
    if (r_stream.eof() || r_stream.bad()) {
        r_stream.clear();
        r_stream.ignore(ignoreSize, '\n');
        return false;
    }

    // Read column index
    r_stream >> i_column;
    if (r_stream.eof() || r_stream.bad()) {
        r_stream.clear();
        r_stream.ignore(ignoreSize, '\n');
        return false;
    }

    // Read value if requested
    r_stream >> value;
    if (r_stream.eof() || r_stream.bad()) {
        r_stream.clear();
        r_stream.ignore(ignoreSize, '\n');
        return false;
    }

    r_stream.ignore(ignoreSize, '\n');

    // Convert (row,column) indices to 0-based indices
    --i_row;
    --i_column;

    // Record entry in the output matrix
    r_matrix.columnIndices.push_back(i_column);
    r_matrix.nonzeros.push_back(value);
    if (r_matrix.rowExtents.size() <= i_row + 1) {
        for (auto i=r_matrix.rowExtents.size(); i<=i_row + 1; ++i) r_matrix.rowExtents.push_back(r_matrix.nonzeros.size());
    }

    return true;
}


TestCSRMatrix parseMatrixMarket(std::istream& r_stream)
{
    TestCSRMatrix matrix = parseMatrixMarketHeader(r_stream);
    while (parseSparseDataLine(r_stream, matrix)) {}
    matrix.rowExtents.push_back(matrix.nonzeros.size());
    return matrix;
}


} // namespace mcgs
