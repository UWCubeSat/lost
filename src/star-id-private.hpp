// Stuff that tests need to access, but not part of the public interface

#ifndef STAR_ID_PRIVATE_H
#define STAR_ID_PRIVATE_H

#include <limits>
#include <utility>
#include <vector>
#include <map>

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


/**
 * Given the result of a pair-distance kvector query, build a hashmultimap of stars to other stars
 * that appeared with it in the query.
 *
 * The resulting map is "symmetrical" in the sense that if a star B is in the map for star A, then
 * star A is also in the map for star B.
 */
std::multimap<int16_t, int16_t> PairDistanceQueryToMap(const int16_t *pairs, const int16_t *end);

/**
 * Try to identify a pattern of starts using pair distances.
 *
 * Pass it 3D spatials, and it will try to match catalog indices
 * @return the "number" of matches. -1 if it cannot possibly match because the inter-star distances are outside the allowable range. 0 if just no match. 1 if unique match (this is the only case you can rely on the values of matchedCatalogIndices). Or, if not a unique match, the number of matches.
 */
template <int numPatternStars>
int IdentifyPatternPairDistance(const PairDistanceKVectorDatabase &db,
                                const Catalog &catalog,
                                float tolerance,
                                const Vec3 spatials[numPatternStars],
                                int matchedCatalogIndices[numPatternStars]) {
    static_assert(numPatternStars >= 3, "Cannot identify patterns w/ less than 3 stars");

    // double check that spatials are unit vectors
    assert(abs(spatials[0].MagnitudeSq() - 1) < 0.0001);

    const Vec3 &iSpatial = spatials[0];
    const Vec3 &jSpatial = spatials[1];
    const Vec3 &kSpatial = spatials[2];

    // sign of determinant, to detect flipped patterns
    bool spectralTorch = iSpatial.CrossProduct(jSpatial)*kSpatial > 0;

    // 2D array of distances between stars
    float distances[numPatternStars][numPatternStars];
    for (int i = 0; i < numPatternStars; i++) {
        for (int j = i+1; j < numPatternStars; j++) {
            distances[i][j] = AngleUnit(spatials[i], spatials[j]);
            distances[j][i] = distances[i][j]; // don't really need symmetry, but why not?

            // the first star is the one we do the actual db queries on, so we want the other distances to be in range
            if (i == 0 &&
                (distances[i][j] < db.MinDistance() + tolerance
                 || distances[i][j] > db.MaxDistance() - tolerance)) {
                return -1;
            }
        }
    }

    std::multimap<int16_t, int16_t> pairDistanceMaps[numPatternStars];
    // we'll build the maps lazily, so record which ones have been built.
    bool builtMaps[numPatternStars] = {false};

    const int16_t *ijEnd, *ikEnd;
    const int16_t *const ijQuery = db.FindPairsExact(catalog, distances[0][1] - tolerance, distances[0][1] + tolerance, &ijEnd);
    const int16_t *const ikQuery = db.FindPairsExact(catalog, distances[0][2] - tolerance, distances[0][2] + tolerance, &ikEnd);
    pairDistanceMaps[2] = PairDistanceQueryToMap(ikQuery, ikEnd);
    builtMaps[2] = true;

    int numMatches = 0;

    Vec3 candidateSpatials[numPatternStars];
    int16_t candidateCatalogIndices[numPatternStars];
    for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd; iCandidateQuery++) {
        candidateCatalogIndices[0] = *iCandidateQuery;
        // depending on parity, the first or second star in the pair is the "other" one
        candidateCatalogIndices[1] = (iCandidateQuery - ijQuery) % 2 == 0
            ? iCandidateQuery[1]
            : iCandidateQuery[-1];

        candidateSpatials[0] = catalog[candidateCatalogIndices[0]].spatial;
        candidateSpatials[1] = catalog[candidateCatalogIndices[1]].spatial;

        Vec3 ijCandidateCross = candidateSpatials[0].CrossProduct(candidateSpatials[1]);

        for (auto kCandidateIt = pairDistanceMaps[2].equal_range(candidateCatalogIndices[0]); kCandidateIt.first != kCandidateIt.second; kCandidateIt.first++) {
            // kCandidate.first is iterator, then ->second is the value (other star)
            candidateCatalogIndices[2] = kCandidateIt.first->second;
            candidateSpatials[2] = catalog[candidateCatalogIndices[2]].spatial;
            bool candidateSpectralTorch = ijCandidateCross*candidateSpatials[2] > 0;
            // checking the spectral-ity early to fail fast
            if (candidateSpectralTorch != spectralTorch) {
#if LOST_DEBUG > 3
                std::cerr << "skipping candidate " << iCandidate << " " << jCandidate << " " << kCandidate << " because spectral-ity mismatch" << std::endl;
#endif
                continue;
            }

            float jkCandidateDist = AngleUnit(candidateSpatials[1], candidateSpatials[2]);
            if (jkCandidateDist < distances[1][2] - tolerance || jkCandidateDist > distances[1][2] + tolerance) {
                continue;
            }

            // now draw the rest of the fucking owl
            for (int l = 3; l < numPatternStars; l++) {
                if (!builtMaps[l]) {
                    const int16_t *lrEnd;
                    const int16_t *const lrQuery = db.FindPairsExact(catalog, distances[l][0] - tolerance, distances[l][0] + tolerance, &lrEnd);
                    pairDistanceMaps[l] = PairDistanceQueryToMap(lrQuery, lrEnd);
                    builtMaps[l] = true;
                }

                auto lCandidateIt = pairDistanceMaps[l].equal_range(candidateCatalogIndices[0]);
                // break fast to the top if there are no matches for this star
                if (lCandidateIt.first == lCandidateIt.second) {
                    // no matches for this star
                    goto nextICandidate;
                }
                for (; lCandidateIt.first != lCandidateIt.second; lCandidateIt.first++) {
                    // lCandidate.first is iterator, then ->second is the value (other star)
                    candidateCatalogIndices[l] = lCandidateIt.first->second;
                    candidateSpatials[l] = catalog[candidateCatalogIndices[l]].spatial;

                    // check distances against all other stars, except 0, because that's part of the query
                    for (int m = 1; m < l; m++) {
                        float lmCandidateDist = AngleUnit(candidateSpatials[l], candidateSpatials[m]);
                        if (lmCandidateDist < distances[l][m] - tolerance || lmCandidateDist > distances[l][m] + tolerance) {
                            goto nextLCandidate;
                        }
                    }

                    // if we get here, we have a match!
                    numMatches++;
                    // if there are multiple matches, the matched indices are undefined anyway, so just set them unconditionally
                    for (int m = 0; m < numPatternStars; m++) {
                        matchedCatalogIndices[m] = candidateCatalogIndices[m];
                    }

                nextLCandidate:;
                }
            }
        }
    nextICandidate:;
    }

    return numMatches;
}

} // namespace lost

#endif
