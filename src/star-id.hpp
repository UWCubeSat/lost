#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "camera.hpp"
#include "centroiders.hpp"
#include "databases.hpp"
#include "star-utils.hpp"

namespace lost {

/**
 * A star idenification algorithm.
 * An algorithm which takes a list of centroids plus some (possibly algorithm-specific) database
 * and then determines which centroids corresponds to which catalog stars.
 */
class StarIdAlgorithm {
   public:
    /// Actualy perform the star idenification. This is the "main" function for StarIdAlgorithm
    virtual StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &,
                               const Camera &) const = 0;

    virtual ~StarIdAlgorithm(){};
};

/**
 * Tetra star identification algorithm based on the one originally proposed by Brown & Stubis
 * Utilizes hashing of 4-star patterns for quick lookup time, decently high accuracy
 */
class TetraStarIdAlgorithm : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &centroids,
                       const Catalog &catalog, const Camera &) const;

   private:
    // Do NOT modify these parameters unless you know what you're doing
    const int numPattStars = 4;
    const int numPattBins = 50;
    const float pattErrorRange = 0.002;
    const float pattMaxError = 0.0005;

};

/// A star-id algorithm that returns random results. For debugging.
class DummyStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &,
                       const Camera &) const;
};

/**
 * A star-id algorithm based on assigning votes for each centroid-catalog pair then choosing the
 * highest voted catalog stars. While the geometric voting algorithm is the simplest true star-id
 * algorithm I know of and is quite fast, it has some reliability issues and there are no
 * statistical guarantees on how often it will return the wrong result. It will also frequently fail
 * to idenify images with few stars.
 */
class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &,
                       const Camera &) const;

    /**
     * @param tolerance Angular tolerance (Two inter-star distances are considered the same if
     * within this many radians)
     */
    explicit GeometricVotingStarIdAlgorithm(float tolerance) : tolerance(tolerance){};

   private:
    float tolerance;
};

/**
 * The "de facto" star-id algorithm used in many real-world missions.
 * Pyramid searches through groups of 4 stars in the image. For each one it tries to find a
 * corresponding 4-star pattern in the database. Four stars is enough that pyramid is often able to
 * uniquely match the first 4-star pattern it tries, making it fast and reliable. However, this only
 * holds true if the camera is calibrated and has low centroid error.
 */
class PyramidStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &,
                       const Camera &) const;
    /**
     * @param tolerance Angular tolerance (Two inter-star distances are considered the same if
     * within this many radians)
     * @param numFalseStars an estimate of the number of false stars in the whole celestial sphere
     * (not just the field of view). Eg, if you estimate 10 dead pixels in a 40 degree FOV, you'd
     * want to multiply that up to a hundred-something numFalseStars.
     * @param maxMismatchProbability The maximum allowable probability for any star to be mis-id'd.
     * @param cutoff Maximum number of pyramids to iterate through before giving up.
     */
    PyramidStarIdAlgorithm(float tolerance, int numFalseStars, float maxMismatchProbability,
                           long cutoff)
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
