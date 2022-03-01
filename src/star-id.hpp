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

class DummyStarIdAlgorithm final : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const override;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const override;
    GeometricVotingStarIdAlgorithm(float tolerance): tolerance(tolerance) { };
private:
    float tolerance;
};


class PyramidStarIdAlgorithm final : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const override;
    /**
     * @param tolerance Angular tolerance in distances (measurement error)
     * @param numFalseStars an estimate of the number of false stars in the whole celestial sphere
     * (not just the field of view). Eg, if you estimate 10 dead pixels in a 40 degree FOV, you'd
     * want to multiply that up to a hundred-something numFalseStars.
     * @param maxMismatchProbability The maximum allowable probability for any star to be mis-id'd.
     * @param cutoff Maximum number of pyramids to iterate through.
     */
    PyramidStarIdAlgorithm(float tolerance, int numFalseStars, float maxMismatchProbability, long cutoff)
        : tolerance(tolerance), numFalseStars(numFalseStars),
          maxMismatchProbability(maxMismatchProbability), cutoff(cutoff) { };
private:
    float tolerance;
    int numFalseStars;
    float maxMismatchProbability;
    long cutoff;
};

class BayesianStarIdAlgorithm final : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const override;

    /**
     * @param tolerance: Angular tolerance in degrees (centroid/measurement error)
     * @param numFalseStars: expected number of false stars throughout the celestial sphere
     * @param softConfidenceThreshold: If the confidence is below this threshold after considering every star in the image, then try to consider negative evidence.
     * @param hardConfidenceThreshold: Do not return the match if below this threshold of confidence.
     * @param admissibleIgnoredProbability: At each step, allow the algorithm to throw out possibilities as long as the probability measure of the discarded possibilities sums to less than admissibleIgnoredProbability
     */
    BayesianStarIdAlgorithm(float tolerance, int numFalseStars,
                            float softConfidenceThreshold, float hardConfidenceThreshold,
                            float admissableIgnoredProbability);

private:
    float tolerance;
    int numFalseStars;
    float softConfidenceThreshold;
    float hardConfidenceThreshold;
    float admissibleIgnoredProbability;
};

}

#endif
