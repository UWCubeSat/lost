#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"
#include "star-utils.hpp"

namespace lost {

class StarIdAlgorithm {
public:
    virtual Stars Go(const void *database, const Stars &) const = 0;
    virtual ~StarIdAlgorithm() { };
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    Stars Go(const void *database, const Stars &) const;
};

class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    Stars Go(const void *database, const Stars &) const;
};

}

#endif
