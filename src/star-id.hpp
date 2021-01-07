#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"
#include "star-utils.hpp"

namespace lost {

class StarIdAlgorithm {
public:
    virtual StarIdentifiers Go(const unsigned char *database, const Stars &) const = 0;
    virtual ~StarIdAlgorithm() { };
};

class DummyStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &) const;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &) const;
    //convert x and y coordinates to degree differences 
    //give a greater range for min-max Query for bigger radius(?)
    //us voting system 
    //optimizations? N^2
    //testing, add false stars and see if the accuracy is still good (maybe just 1 or 2 false stars)
    //instead of random testing, maybe we have random stars from our db then convert into pixel values
};


class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &) const;
};

}

#endif
