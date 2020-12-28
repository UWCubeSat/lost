#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"
#include "star-utils.hpp"

namespace lost {

class StarIdAlgorithm {
public:
    virtual void Go(const unsigned char *database, Stars *) const = 0;
    virtual ~StarIdAlgorithm() { };
};

class DummyStarIdAlgorithm : public StarIdAlgorithm {
public:
    void Go(const unsigned char *database, Stars *) const;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    void Go(const unsigned char *database, Stars *) const;
};

class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    void Go(const unsigned char *database, Stars *) const;
};

}

#endif
