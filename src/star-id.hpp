#ifndef STAR_ID_H
#define STAR_ID_H

#define _USE_MATH_DEFINES

#include <cmath>  // TODO: added in here
#include <vector>

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "centroiders.hpp"
#include "star-utils.hpp"

using namespace std;

namespace lost {


class StarIdAlgorithm {
   public:
    virtual StarIdentifiers Go(const unsigned char *database, const Stars &,
                               const Catalog &, const Camera &) const = 0;
    virtual ~StarIdAlgorithm(){};
};


class DummyStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
    GeometricVotingStarIdAlgorithm(float tolerance) : tolerance(tolerance){};

   private:
    float tolerance;
};

class PyramidStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
    /**
     * @param tolerance Angular tolerance in distances (measurement error)
     * @param numFalseStars an estimate of the number of false stars in the
     * whole celestial sphere (not just the field of view). Eg, if you estimate
     * 10 dead pixels in a 40 degree FOV, you'd want to multiply that up to a
     * hundred-something numFalseStars.
     * @param maxMismatchProbability The maximum allowable probability for any
     * star to be mis-id'd.
     * @param cutoff Maximum number of pyramids to iterate through.
     */
    PyramidStarIdAlgorithm(float tolerance, int numFalseStars,
                           float maxMismatchProbability, long cutoff)
        : tolerance(tolerance),
          numFalseStars(numFalseStars),
          maxMismatchProbability(maxMismatchProbability),
          cutoff(cutoff){};

   private:
    float tolerance;
    int numFalseStars;
    float maxMismatchProbability;
    long cutoff;
};

}  // namespace lost

#endif
