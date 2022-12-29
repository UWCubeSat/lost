#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "camera.hpp"
#include "centroiders.hpp"
#include "star-utils.hpp"
#include "databases.hpp"

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

class TetraStarIdAlgorithm : public StarIdAlgorithm {
 public:
  StarIdentifiers Go(const unsigned char *database, const Stars &centroids, const Catalog &catalog,
                     const Camera &) const;

 private:
  // const float fov = 25.5705;   // in degrees, TODO: not used anywhere so delete later
  // TODO: this should be read from the database in the Go algorithm
  // TODO: update database for this
  const float maxFov = 12.00;  // in degrees, max FOV of database
  // TODO: this may not be accurate, think I saw a 20 FOV somewhere. Also make this part of
  // constructor / default, not hardcoded

  // I feel these should be held constant, cannot be changed
  const int numPattStars = 4;
  const int numPattBins = 25;
  const float pattMaxError = 0.005;

  // const int catalogLength = 11841082; // default database
  // const int catalogLength = 8979154; // tetra3 fov=12, stable
  // number of patterns in catalog:
  // const int catalogLength = 8978892;  // hardcoded, just for testing - remove later
  const int catalogLength = 8951660;

  const long long MAGIC_RAND = 2654435761;

  /**
   * @brief Hash function, convert from star pattern representation into index in pattern catalog
   *
   * A pattern is represented by a sorted, 5-element vector of edge ratios (divided by largest edge)
   * Technically of dimension C(numPattStars, 2) - 1
   *
   * @param key 5-element vector, sorted, of pattern edge ratios
   * @param binFactor
   * @param maxIndex Number of rows in pattern catalog
   * @return int
   */
  int KeyToIndex(std::vector<int> key, int binFactor, int maxIndex) const;

  /**
   * @brief Get all possible matching patterns starting from given index
   *
   * Perform quadratic probing
   *
   * @param index
   * @param pattCatFile
   * @return std::vector<std::vector<int>> List of 4-star patterns that could be matches
   */
  // std::vector<std::vector<int>> GetAtIndex(int index, std::ifstream &pattCatFile) const;
    std::vector<std::vector<int>> GetAtIndex(int index, int maxIndex, const TetraDatabase &db) const;
  // TODO: change, should read from database not the file

  // std::vector<std::vector<int>> GetAtIndex(int index, TetraDatabase db) const; REMOVE
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
   * @param tolerance Angular tolerance (Two inter-star distances are considered the same if within
   * this many radians)
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
   * @param tolerance Angular tolerance (Two inter-star distances are considered the same if within
   * this many radians)
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
