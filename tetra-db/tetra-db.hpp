#ifndef TETRA_DB_H
#define TETRA_DB_H

#define _USE_MATH_DEFINES
#include "math.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <map>
#include <array>
#include <set>

struct Bsc5DataEntry {
    double ra_1950;
    double dec_1950;
    float id;
    char type[2];
    short magnitude;
    float ra_pm;
    float dec_pm;
};

struct StarTableEntry {
    float ra;
    float dec;
    float x;
    float y;
    float z;
    float magnitude;
    float id; // really int
};

typedef std::array<short, 3> ShortVec3;
typedef std::array<float, 3> FloatVec3;

class Tetra {
   public:
    void static GenerateDatabase(float maxFOV, std::string &starTableFName,
                          std::string &catalogFName, std::string &propsFName,
                          int pattStarsPerFOV = 10,
                          int verificationStarsPerFOV = 20,
                          float starMaxMagnitude = 7,
                          float starMinSeparation = 0.05,
                          float pattMaxError = 0.005);

    int static KeyToIndex(std::vector<int> key, int binFactor, int maxIndex);
};

#endif
