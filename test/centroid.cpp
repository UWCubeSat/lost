#include <catch.hpp>

#include <math.h>
#include <iostream>

namespace lost {
    int RowOrColumn(long i, int leftover, int div);
    int FindSubdivision(long i, int imageWidth, int imageHeight, int subdivisions);
    int StartOfSubdivision(int i, int leftover, int div);
}


TEST_CASE("Testing Start and End Rows") {
    int row = 0;
    int horizontalLeftover = 5;
    int horizontalDiv = 10;
    CHECK(lost::StartOfSubdivision(row, horizontalLeftover, horizontalDiv) == 0);
    CHECK(lost::StartOfSubdivision(row + 1, horizontalLeftover, horizontalDiv) == 11);
}

TEST_CASE("Testing Start and End Columns") {
    int col = 0;
    int verticalLeftover = 5;
    int verticalDiv = 10;
    CHECK(lost::StartOfSubdivision(col, verticalLeftover, verticalDiv) == 0);
    CHECK(lost::StartOfSubdivision(col + 1, verticalLeftover, verticalDiv) == 11);
}

TEST_CASE("Correct Row: 4th Row") {
    CHECK(lost::RowOrColumn(33795 / 1024, 1024 % 100, 1024 / 100) == 3);
}

TEST_CASE("Correct Row: 27th Row") {
    CHECK(lost::RowOrColumn((24 * 1024 * 11 + 10 * 2 * 1024 + 20) / 1024, 1024 % 100, 1024 / 100) == 26);
}

TEST_CASE("Correct Row: Last Row 1") {
    CHECK(lost::RowOrColumn((1024 * 1024 - 1) / 1024, 1024 % 100, 1024 / 100) == 99);
}

TEST_CASE("Correct Box: Last Subdivision") {
    CHECK(lost::FindSubdivision(1024 * 1024 - 1, 1024, 1024, 1024) == 1024 * 1024 - 1);
}

TEST_CASE("Testing Inverse Conjecture") {
    int size = 1024;
    int subdivisions = 10;
    int leftover = size % subdivisions;
    int div = size / subdivisions;
    for(int i = 0; i < subdivisions; i++) {
        CHECK(lost::RowOrColumn(lost::StartOfSubdivision(i, leftover, div), leftover, div) == i);
    }
}

