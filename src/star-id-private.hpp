// Stuff that tests need to access, but not part of the public interface

#ifndef STAR_ID_PRIVATE_H
#define STAR_ID_PRIVATE_H

#include "star-id.hpp"
#include "databases.hpp"

#include <limits>

namespace lost {

/// unidentified centroid used in IdentifyRemainingStarsPairDistance
/// The "angles" through here are "triangular angles". A triangular angle is a 2D angle in a triangle formed by three centroids.
class IRUnidentifiedCentroid {
public:
    IRUnidentifiedCentroid(const Star &star)
        : bestAngleFrom90(std::numeric_limits<float>::max()), // should be infinity
          bestStar1(0,0), bestStar2(0,0),
          star(&star) { }

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
                                       const Camera &,
                                       float tolerance);

}

#endif
