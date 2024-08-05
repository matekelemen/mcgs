#pragma once

// --- STL Includes ---
#include <vector> // std::vector
#include <istream> // std::istream
#include <variant> // std::variant


namespace mcgs {


struct TestCSRMatrix
{
    using Index = unsigned long;

    using Value = double;

    Index rowCount;

    Index columnCount;

    Index entryCount;

    std::vector<Index> rowExtents;

    std::vector<Index> columnIndices;

    std::vector<Value> entries;
}; // struct CSRMatrix


using TestDenseVector = std::vector<double>;


std::variant<
    TestCSRMatrix,
    TestDenseVector
> parseMatrixMarket(std::istream& r_stream);


} // namespace mcgs
