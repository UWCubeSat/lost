// Stuff that tests need to access, but not part of the public interface

#ifndef STAR_ID_PRIVATE_H
#define STAR_ID_PRIVATE_H

#include <limits>
#include <utility>
#include <vector>

#include "star-utils.hpp"
#include "star-id.hpp"
#include "databases.hpp"

namespace lost {

/// unidentified centroid used in IdentifyRemainingStarsPairDistance
/// The "angles" through here are "triangular angles". A triangular angle is a 2D angle in a triangle formed by three centroids.
class IRUnidentifiedCentroid {
public:
    IRUnidentifiedCentroid(const Star &star, int16_t index)
        : bestAngleFrom90(std::numeric_limits<float>::max()), // should be infinity
          bestStar1(0,0), bestStar2(0,0),
          index(index),
          star(&star) {
        identifiedStarsInRange.reserve(10); // this does quite measurably improve performance, at least on desktop
    }

    // "null" has index=-1
    IRUnidentifiedCentroid()
        : bestStar1(0,0), bestStar2(0,0),
          index(-1) { }

    float bestAngleFrom90; /// For the pair of other centroids forming the triangular angle closest to 90 degrees, how far from 90 degrees it is (in radians)
    StarIdentifier bestStar1; /// One star corresponding to bestAngleFrom90
    StarIdentifier bestStar2; /// The other star corresponding to bestAngleFrom90
    int16_t index; /// Index into list of all centroids
    const Star *star;

private:
    // possible improvement: Use a tree map here to allow binary search
    std::vector<std::pair<float, StarIdentifier>> identifiedStarsInRange;

private:
    float VerticalAnglesToAngleFrom90(float v1, float v2);

public:
    void AddIdentifiedStar(const StarIdentifier &starId, const Stars &stars);
};

std::vector<int16_t> IdentifyThirdStar(const PairDistanceKVectorDatabase &db,
                                       const Catalog &catalog,
                                       int16_t catalogIndex1, int16_t catalogIndex2,
                                       float distance1, float distance2,
                                       float tolerance);

int IdentifyRemainingStarsPairDistance(StarIdentifiers *,
                                       const Stars &,
                                       const PairDistanceKVectorDatabase &,
                                       const Catalog &,
                                       const Camera &,
                                       float tolerance);

// PYRAMID

/**
 * The best pyramid "starting from" a certain star. The "other" stars are ordered by their distance
 * from the main star for this struct.
 *
 * "Distance" in this class is angular distance, as usual.
 *
 * If distancesSum is nonpositive, no suitable pyramid exists for this star.
 */
class BestPyramidAtStar {
public:
    int16_t centroidIndices[4];

    float distancesSum;

    BestPyramidAtStar(int16_t centroidIndex0, int16_t centroidIndex1, int16_t centroidIndex2, int16_t centroidIndex3,
                      float distancesSum)
        : distancesSum(distancesSum) {
        centroidIndices[0] = centroidIndex0;
        centroidIndices[1] = centroidIndex1;
        centroidIndices[2] = centroidIndex2;
        centroidIndices[3] = centroidIndex3;
    }

    // "no suitable pyramid" constructor
    explicit BestPyramidAtStar(int16_t mainCentroidIndex)
        : distancesSum(-1) {
        centroidIndices[0] = mainCentroidIndex;
        centroidIndices[1] = -1;
        centroidIndices[2] = -1;
        centroidIndices[3] = -1;
    }

    bool operator<(const BestPyramidAtStar &other) const {
        return distancesSum > 0 && distancesSum < other.distancesSum;
    }

    bool isNull() const {
        return distancesSum <= 0;
    }
};

/**
 * Keep finding the best pyramid to attempt to identify next.
 *
 * Rough strategy is to find the overall best pyramid first, then remove all the stars in that
 * pyramid and find the next best one (the idea being that if identification failed for the first
 * pyramid, there is likely a false star in it, so we want to try and identify different stars next
 * for highest chance of success).
 *
 * Of course, it's possible that we'll exhaust all stars using this strategy and still won't have
 * found a valid pyramid, even though valid pyramids totally exist. In that case, we'll just resort
 * to doing something worse TODO.
 */
class PyramidIterator {
public:
    /// Please ensure that `centroids` outlives the PyramidIterator! Also, minDistance and maxDistance are exact, we don't offset by tolerance, you probably want to do that!
    PyramidIterator(const std::vector<Vec3> &centroidSpatials, float minDistance, float maxDistance);

    /// Returns the next best pyramid, or a "no pyramid" pyramid.
    BestPyramidAtStar Next();

private:
    float minDistance, maxDistance;
    // once length of this is less than 4, we switch to alternate strategy:
    const std::vector<Vec3> &allCentroidSpatials;
    std::vector<int16_t> untriedCentroidIndices;
};

enum class IdentifyPyramidResult {
    /// Won't even attempt identification of the given pyramid, not possible to identify for some reason (eg inter-star distances out of k-vector range)
    CannotPossiblyIdentify = -1,
    NoMatch = 0, /// Did not find any match
    MatchedUniquely = 1, /// Matched uniquely, this is the only "successful" result
    MatchedAmbiguously = 2, /// Multiple matches, you shouldn't use any of them
};

/**
 * Try to identify a pattern of four stars.
 *
 * Pass it four 3D spatials, and it will try to find 4 catalog stars which match.
 * @param a,b,c,d Catalog indices of the identified stars, if uniquely matched.
 * @return see IdentifyPyramidResult. You can also sorta treat it like a number, though: "how many matches" were there?
 */
IdentifyPyramidResult IdentifyPyramid(const PairDistanceKVectorDatabase &,
                                      const Catalog &,
                                      float tolerance,
                                      const Vec3 &, const Vec3 &, const Vec3 &, const Vec3 &,
                                      int *a, int *b, int *c, int *d);

} // namespace lost

#endif
