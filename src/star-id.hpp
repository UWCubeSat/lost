#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"

namespace lost {

class IdentifiedStar {
    
};

typedef std::vector<IdentifiedStar> Stars;

class StarIdAlgorithm {
public:
    virtual Stars Go(const void *database, const Centroids &) const = 0;
    virtual ~StarIdAlgorithm() { };
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    Stars Go(const void *database, const Centroids &) const;
};

class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    Stars Go(const void *database, const Centroids &) const;
};

}

#endif
