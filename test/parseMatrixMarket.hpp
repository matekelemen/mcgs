#pragma once

// --- STL Includes ---
#include <vector> // std::vector
#include <fstream> // std::ifstream


namespace mcgs {


struct TestCSRMatrix
{
    using Index = unsigned long;

    using Value = double;

    Index rowCount;

    Index columnCount;

    Index nonzeroCount;

    std::vector<Index> rowExtents;

    std::vector<Index> columnIndices;

    std::vector<Value> nonzeros;
}; // struct CSRMatrix


TestCSRMatrix parseMatrixMarket(std::istream& r_stream);


} // namespace mcgs
