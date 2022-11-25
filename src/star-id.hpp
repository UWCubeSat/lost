#ifndef STAR_ID_H
#define STAR_ID_H

#include <vector>

#include "centroiders.hpp"
#include "star-utils.hpp"
#include "camera.hpp"

namespace lost {

/**
 * A star idenification algorithm.
 * An algorithm which takes a list of centroids plus some (possibly algorithm-specific) database, and then determines which centroids corresponds to which catalog stars.
 */
class StarIdAlgorithm {
public:
    /// Actualy perform the star idenification. This is the "main" function for StarIdAlgorithm
    virtual StarIdentifiers Go(
        const unsigned char *database, const Stars &, const Catalog &, const Camera &) const = 0;

    virtual ~StarIdAlgorithm() { };
};

/// Class just to hold the Tetra "databases" for now
// TODO: remove
// class TetraDatabase {
//    public:
//     // Generated at FOV = 12
//     TetraDatabase()
//         : pattCatalog(11841082, std::vector<int>(4)),
//           starTable(8416, std::vector<float>(7)) {}

//     // pattCatalog is a 2D matrix, with 11841082 rows and 4 columns (starIDs)
//     std::vector<std::vector<int>> pattCatalog;
//     // std::vector<std::vector<int>> pattCatalog (std::vector<int>(4),
//     // 11841082); int pattCatalog[11841082][4];
//     // TODO: seg faulting if I make the vector back into an array
//     // This overflows the stack?

//     // starTable is a 2D matrix, with 8416 rows and 7 columns
//     // Columns: RA, DE, x, y, z, Magnitude, Star ID
//     // x, y, and z correspond to the Vec3 spatial vector for this star
//     // Star ID in each row is exactly what is displayed on annotated.png
//     // TODO: HR number I think?
//     std::vector<std::vector<float>> starTable;
//     // float starTable[8416][7];
//     // TODO: do it this way, otherwise we get a seg fault

//     void fillPattCatalog();
//     void fillStarTable();
// };

class TetraStarIdAlgorithm: public StarIdAlgorithm{
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;

private:
    const float fov = 25.5705; // in degrees
    const float maxFov = 12.00; // in degrees, max FOV of database
    // TODO: this may not be accurate, think I saw a 20 FOV somewhere. Also make this part of constructor / default, not hardcoded

    // I feel these should be held constant, cannot be changed
    const int numPattStars = 4;
    const int numPattBins = 25;
    const float pattMaxError = 0.005;

    // const int catalogLength = 11841082; // default database
    // const int catalogLength = 8979154; // tetra3 fov=12, stable
    // number of patterns in catalog:
    const int catalogLength = 8978892; // hardcoded, just for testing - remove later

    const long long MAGIC_RAND = 2654435761;

    const int starTableRowSize = 7; // TODO: remove

    int KeyToIndex(std::vector<int> key, int binFactor, int maxIndex) const;
    // std::vector<std::vector<int>> GetAtIndex(int index, TetraDatabase db) const;
    std::vector<std::vector<int>> GetAtIndex(int index, std::ifstream &pattCatFile) const;
};

/// A star-id algorithm that returns random results. For debugging.
class DummyStarIdAlgorithm final : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
};

/**
 * A star-id algorithm based on assigning votes for each centroid-catalog pair then choosing the highest voted catalog stars.
 * While the geometric voting algorithm is the simplest true star-id algorithm I know of and is quite fast, it has some reliability issues and there are no statistical guarantees on how often it will return the wrong result. It will also frequently fail to idenify images with few stars.
 */
class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;

    /**
     * @param tolerance Angular tolerance (Two inter-star distances are considered the same if within this many radians)
     */
    explicit GeometricVotingStarIdAlgorithm(float tolerance): tolerance(tolerance) { };
private:
    float tolerance;
};


/**
 * The "de facto" star-id algorithm used in many real-world missions.
 * Pyramid searches through groups of 4 stars in the image. For each one it tries to find a corresponding 4-star pattern in the database. Four stars is enough that pyramid is often able to uniquely match the first 4-star pattern it tries, making it fast and reliable. However, this only holds true if the camera is calibrated and has low centroid error.
 */
class PyramidStarIdAlgorithm final : public StarIdAlgorithm {
public:
    StarIdentifiers Go(const unsigned char *database, const Stars &, const Catalog &, const Camera &) const;
    /**
     * @param tolerance Angular tolerance (Two inter-star distances are considered the same if within this many radians)
     * @param numFalseStars an estimate of the number of false stars in the whole celestial sphere
     * (not just the field of view). Eg, if you estimate 10 dead pixels in a 40 degree FOV, you'd
     * want to multiply that up to a hundred-something numFalseStars.
     * @param maxMismatchProbability The maximum allowable probability for any star to be mis-id'd.
     * @param cutoff Maximum number of pyramids to iterate through before giving up.
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

}

#endif
