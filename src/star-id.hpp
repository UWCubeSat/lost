#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"
#include "star-utils.hpp"
#include "camera.hpp"

namespace lost {

class StarIdAlgorithm {
public:
    virtual StarIdentifiers Go(
        const unsigned char *database, const Stars &, const Catalog &, const Camera &) const = 0;
    virtual ~StarIdAlgorithm() { };
};

class DummyStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
    GeometricVotingStarIdAlgorithm(const float tolerance): tolerance(tolerance) { };
private:
    float tolerance;
};


class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
};

class NonDimStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
    NonDimStarIdAlgorithm(const float tolerance): tolerance(tolerance) { };
private:
    float tolerance;
};

}

#endif
