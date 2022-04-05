#include <catch.hpp>

#include <math.h>
#include <iostream>


int start(int imageWidth, int imageHeight, int subdivisions, int i) {
    int leftover = imageHeight % subdivisions;
    int div = imageHeight / subdivisions;
    // int i = 1;
    // int start = 0;
    // int max = 0;
    if(i < leftover) {
        return i * (div + 1) * imageWidth;
        // max = (i+1) * (div + 1) * imageWidth;
    } else {
        return leftover * (div + 1) * imageWidth + (i - leftover) * (div) * imageWidth;
        // max = leftover * (div + 1) * imageWidth + (i + 1 - leftover) * (div) * imageWidth;
    }
}


TEST_CASE("Sanity Check") {
    CHECK(start(1024, 1024, 100, 0) == 0);
}

TEST_CASE("Start: 3rd Element") {
    CHECK(start(1024,1024,100, 2) == 22*1024);
}

TEST_CASE("Start: 27th Element") {
    CHECK(start(1024,1024,100, 26) == (11*24+10*2)*1024);
}
