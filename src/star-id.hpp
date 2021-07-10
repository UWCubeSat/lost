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
    GeometricVotingStarIdAlgorithm(float tolerance): tolerance(tolerance) { };
private:
    float tolerance;
};


class PyramidStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
    // tolerance is an angular distances for two angles to be considered the same, maxProbability is
    // the maximum allowable likelihood for a match to have been found in a uniformly random star
    // layout, and cutoff is the maximum number of pyramids to inspect before giving up (if there
    // are dozens of stars in the field of view, inspecting all pyramids will be very slow)
    PyramidStarIdAlgorithm(float tolerance, float maxProbability, long cutoff)
        : tolerance(tolerance), maxProbability(maxProbability), cutoff(cutoff) { };
private:
    float tolerance;
    float maxProbability;
    long cutoff;
};

}

#endif
