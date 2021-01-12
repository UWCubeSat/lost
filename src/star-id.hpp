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
};


class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &) const;
};

}

#endif
