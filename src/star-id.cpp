#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <utility>
#include <map>

#include "star-id.hpp"
#include "star-id-private.hpp"
#include "databases.hpp"
#include "attitude-utils.hpp"

namespace lost {

StarIdentifiers DummyStarIdAlgorithm::Go(
    const unsigned char *, const Stars &stars, const Catalog &catalog, const Camera &) const {

    StarIdentifiers result;

    unsigned int randomSeed = 123456;
    for (int i = 0; i < (int)stars.size(); i++) {
        result.push_back(StarIdentifier(i, rand_r(&randomSeed) % catalog.size()));
    }

    return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    for (int i = 0; i < (int)stars.size(); i++) {
        std::vector<int16_t> votes(catalog.size(), 0);
        Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                // TODO: find a faster way to do this:
                std::vector<bool> votedInPair(catalog.size(), false);
                Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
                float greatCircleDistance = AngleUnit(iSpatial, jSpatial);
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = greatCircleDistance - tolerance;
                float upperBoundRange = greatCircleDistance + tolerance;
                const int16_t *upperBoundSearch;
                const int16_t *lowerBoundSearch = vectorDatabase.FindPairsLiberal(
                    lowerBoundRange, upperBoundRange, &upperBoundSearch);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastructure
                for (const int16_t *k = lowerBoundSearch; k != upperBoundSearch; k++) {
                    if ((k - lowerBoundSearch) % 2 == 0) {
                        float actualAngle = AngleUnit(catalog[*k].spatial, catalog[*(k+1)].spatial);
                        assert(actualAngle <= greatCircleDistance + tolerance * 2);
                        assert(actualAngle >= greatCircleDistance - tolerance * 2);
                    }
                    if (!votedInPair[*k] || true) {
                        // if (i == 542 && *k == 9085) {
                        //     printf("INC, distance %f from query %f to %f\n", greatCircleDistance,
                        //         lowerBoundRange, upperBoundRange);
                        // }
                        votes[*k]++;
                        votedInPair[*k] = true;
                    }
                }
                // US voting system
            }
        }
        // Find star w most votes
        int16_t maxVotes = votes[0];
        int indexOfMax = 0;
        for (int v = 1; v < (int)votes.size(); v++) {
            if (votes[v] > maxVotes) {
                maxVotes = votes[v];
                indexOfMax = v;
            }
        }
        // if (i == 542) {
        //     for (float dist : vectorDatabase.StarDistances(9085, catalog)) {
        //         printf("Actual 9085 distance: %f\n", dist);
        //     }
        //     puts("Debug star.");
        //     for (int i = 0; i < (int)votes.size(); i++) {
        //         if (votes[i] > maxVotes/2) {
        //             printf("Star %4d received %d votes.\n", catalog[i].name, votes[i]);
        //         }
        //     }
        //     printf("Debug star: Actually voted for %d with %d votes\n",
        //            catalog[indexOfMax].name, maxVotes);
        // }
        // printf("Max votes: %d\n", maxVotes);
        //starIndex = i, catalog index = indexOfMax
        StarIdentifier newStar(i, indexOfMax);
        // Set identified[i] to value of catalog index of star w most votesr
        identified.push_back(newStar);
    }
    //optimizations? N^2
    //https://www.researchgate.net/publication/3007679_Geometric_voting_algorithm_for_star_trackers
    //
    // Do we have a metric for localization uncertainty? Star brighntess?
    //loop i from 1 through n
    std::vector<int16_t> verificationVotes(identified.size(), 0);
    for (int i = 0; i < (int)identified.size(); i++) {
        //loop j from i+1 through n
        for (int j = i + 1; j < (int)identified.size(); j++) {
            // Calculate distance between catalog stars
            CatalogStar first = catalog[identified[i].catalogIndex];
            CatalogStar second = catalog[identified[j].catalogIndex];
            float cDist = AngleUnit(first.spatial, second.spatial);

            Star firstIdentified = stars[identified[i].starIndex];
            Star secondIdentified = stars[identified[j].starIndex];
            Vec3 firstSpatial = camera.CameraToSpatial(firstIdentified.position);
            Vec3 secondSpatial = camera.CameraToSpatial(secondIdentified.position);
            float sDist = Angle(firstSpatial, secondSpatial);

            //if sDist is in the range of (distance between stars in the image +- R)
            //add a vote for the match
            if (abs(sDist - cDist) < tolerance) {
                verificationVotes[i]++;
                verificationVotes[j]++;
            }
        }
    }
    // Find star w most votes
    int maxVotes = verificationVotes.size() > 0 ? verificationVotes[0] : 0;
    for (int v = 1; v < (int)verificationVotes.size(); v++) {
        if (verificationVotes[v] > maxVotes) {
            maxVotes = verificationVotes[v];
        }
    }

    // If the stars are within a certain range of the maximal number of votes,
    // we consider it correct.
    // maximal votes = maxVotes
    StarIdentifiers verified;
    int thresholdVotes = maxVotes * 3 / 4;
    printf("Verification threshold: %d\n", thresholdVotes);
    for (int i = 0; i < (int)verificationVotes.size(); i++) {
        if (verificationVotes[i] > thresholdVotes) {
            verified.push_back(identified[i]);
        }
    }

    return verified;
}

/*
 * Strategies:
 * 1. For each star, enumerate all stars which have the same combination of distances to some
 *  other stars, getting down to a hopefully small (<10) list of candidates for each star, then
 *  do a quad-nested loop to correlate them.
 * 2. Loop through all possible stars in the catalog for star i. Then look at edge ij, using
 * this to select possible j-th stars. If ever there is not a possible j-th star, continue the
 * i-loop. When a possible ij combination is found, loop through k stars according to ik. IF
 * none are found, continue the outer i loop. If some are found, check jk for each one. For each possible ijk triangle,
 */

/**
 * Given a list of star pairs, finds all those pairs which involve a certain star.
 * Here "involve" means that one of the two stars in the pair is the given star.
 */
class PairDistanceInvolvingIterator {
public:
    /**
     * Create a "past-the-end" iterator.
     * If another PairDistanceInvolvingIterator is equal to this, then it is done iterating.
     */
    PairDistanceInvolvingIterator()
        : pairs(NULL), end(NULL) { };

    /**
     * The main constructor.
     * @param pairs Start of an array of star pairs
     * @param end "past-the-end" pointer for \p pairs
     * @param The catalog index that we want to be involved in the outputted pairs.
     */
    PairDistanceInvolvingIterator(const int16_t *pairs, const int16_t *end, int16_t involving)
        : pairs(pairs), end(end), involving(involving) {

        assert((end-pairs)%2 == 0);
        ForwardUntilInvolving();
    };

    // PairDistanceInvolvingIterator operator++() {
    //     PairDistanceInvolvingIterator result(*this);
    //     ++(*this);
    //     return result;
    // }

    /// Move to the next matching pair.
    PairDistanceInvolvingIterator &operator++() {
        assert(HasValue());
        pairs += 2;
        ForwardUntilInvolving();
        return *this;
    }

    /// Access the curent pair.
    int16_t operator*() const {
        return curValue;
    }

    /// Whether the iterator is currently on a value. (false if iteration is complete)
    bool HasValue() {
        return pairs != end;
    }

    // bool operator==(const PairDistanceInvolvingIterator &other) const {
    //     return ()other.pairs == pairs;
    // }

    // bool operator!=(const PairDistanceInvolvingIterator &other) const {
    //     return !(*this == other);
    // }
private:
    const int16_t *pairs;
    const int16_t *end;
    int16_t involving;
    int16_t curValue;

    /// like postfix++, except it's a no-op if already on a valid spot.
    void ForwardUntilInvolving() {
        while (pairs != end) {
            if (pairs[0] == involving) {
                curValue = pairs[1];
                return;
            }
            if (pairs[1] == involving) {
                curValue = pairs[0];
                return;
            }
            pairs += 2;
        }
    }
};

std::vector<int16_t> ConsumeInvolvingIterator(PairDistanceInvolvingIterator it) {
    std::vector<int16_t> result;
    for (; it.HasValue(); ++it) {
        result.push_back(*it);
    }
    return result;
}

/**
 * Given the result of a pair-distance kvector query, build a hashmultimap of stars to other stars
 * that appeared with it in the query.
 * 
 * The resulting map is "symmetrical" in the sense that if a star B is in the map for star A, then
 * star A is also in the map for star B.
 */
std::multimap<int16_t, int16_t> PairDistanceQueryToMap(const int16_t *pairs, const int16_t *end) {
#if LOST_DEBUG > 3
    std::cerr << "Building map from " << (end-pairs)/2 << " pairs" << std::endl;
#endif
    std::multimap<int16_t, int16_t> result;
    // result.reserve(end - pairs); // equals number of pairs times two, which is correct.
    for (const int16_t *p = pairs; p != end; p += 2) {
        result.emplace(p[0], p[1]);
        result.emplace(p[1], p[0]);
    }
    return result;
}

float IRUnidentifiedCentroid::VerticalAnglesToAngleFrom90(float v1, float v2) {
    return abs(FloatModulo(v1-v2, M_PI) - M_PI_2);
}

/**
 * When a centroid within range of this centroid is identified, call this function. This function
 * does /not/ check whether the centroid is within range.
 */
void IRUnidentifiedCentroid::AddIdentifiedStar(const StarIdentifier &starId, const Stars &stars) {
    const Star &otherStar = stars[starId.starIndex];
    Vec2 positionDifference = otherStar.position - star->position;
    float angleFromVertical = atan2(positionDifference.y, positionDifference.x);

    for (const auto &otherPair : identifiedStarsInRange) {
        float curAngleFrom90 = VerticalAnglesToAngleFrom90(otherPair.first, angleFromVertical);
        if (curAngleFrom90 < bestAngleFrom90) {
            bestAngleFrom90 = curAngleFrom90;
            bestStar1 = starId;
            bestStar2 = otherPair.second;
        }
    }

    identifiedStarsInRange.emplace_back(angleFromVertical, starId);
}

/**
 * Return all the unidentified centroids within the requested distance bounds from `star`
 *
 * The returned vector has pointers into the vector passed as an argument. Thus, it's important not
 * to modify the `centroids` argument after calling.
 */
std::vector<std::vector<IRUnidentifiedCentroid *>::iterator> FindUnidentifiedCentroidsInRange(
    std::vector<IRUnidentifiedCentroid *> *centroids, const Star &star, const Camera &camera,
    float minDistance, float maxDistance) {

    Vec3 ourSpatial = camera.CameraToSpatial(star.position).Normalize();

    float minCos = cos(maxDistance);
    float maxCos = cos(minDistance);

    std::vector<std::vector<IRUnidentifiedCentroid *>::iterator> result;
    for (auto it = centroids->begin(); it != centroids->end(); ++it) {
        Vec3 theirSpatial = camera.CameraToSpatial((*it)->star->position).Normalize();
        float angleCos = ourSpatial * theirSpatial;
        if (angleCos >= minCos && angleCos <= maxCos) {
            result.push_back(it);
        }
    }

    return result;

    // TODO: optimize by sorting on x-coordinate (like in tracking), or maybe even kd-tree

    // // Find the first centroid that is within range of the given centroid.
    // auto firstInRange = std::lower_bound(centroids->begin(), centroids->end(), star.position.x - maxDistance,
    //     [](const IRUnidentifiedCentroid &centroid, float x) {
    //         return centroid.star.position.x < x;
    //     });

    // // Find the first centroid that is not within range of the given centroid.
    // auto firstNotInRange = std::lower_bound(firstInRange, centroids->end(), star.position.x + maxDistance,
    //     [](const IRUnidentifiedCentroid &centroid, float x) {
    //         return centroid.star.position.x <= x;
    //     });

    // // Copy the pointers to the stars into the result vector.
    // for (auto it = firstInRange; it != firstNotInRange; ++it) {
    //     float distance = Distance(star.position, it->star.position);
    //     if (distance >= minDistance && distance <= maxDistance) {
    //         result.push_back(&*it);
    //     }
    // }
}

/**
 * Given a list of unidentified centroids not yet at the soft threshold, and a list of unidentified
 * centroids already below the soft threshold, appropriately add the given centroid to all the
 * unidentified centroids still above the threshold, and perhaps move them to the below threshold
 * list.
 *
 * @param angleFrom90Threshold Once an IRUnidentifiedCentroid's best angle from 90 goes below this threshold
 */
void AddToAllUnidentifiedCentroids(const StarIdentifier &starId, const Stars &stars,
                                   std::vector<IRUnidentifiedCentroid *> *aboveThresholdCentroids,
                                   std::vector<IRUnidentifiedCentroid *> *belowThresholdCentroids,
                                   float minDistance, float maxDistance,
                                   float angleFrom90Threshold,
                                   const Camera &camera) {

    std::vector<int16_t> nowBelowThreshold; // centroid indices newly moved above the threshold
    // don't need to iterate through the centroids that are already below the threshold, for performance.
    for (auto centroidIt : FindUnidentifiedCentroidsInRange(aboveThresholdCentroids, stars[starId.starIndex], camera, minDistance, maxDistance)) {
        (*centroidIt)->AddIdentifiedStar(starId, stars);
        if ((*centroidIt)->bestAngleFrom90 <= angleFrom90Threshold) {
            belowThresholdCentroids->push_back(*centroidIt);
            nowBelowThreshold.push_back((*centroidIt)->index);
        }
    }
    // remove all centroids with indices in nowBelowThreshold from aboveThresholdCentroids
    aboveThresholdCentroids->erase(std::remove_if(aboveThresholdCentroids->begin(), aboveThresholdCentroids->end(),
        [&nowBelowThreshold](const IRUnidentifiedCentroid *centroid) {
            return std::find(nowBelowThreshold.begin(), nowBelowThreshold.end(), centroid->index) != nowBelowThreshold.end();
        }), aboveThresholdCentroids->end());
}

/**
 * Given two already-identified centroids, and the distance from each to an as-yet unidentified
 * third centroid, return a list of candidate catalog stars that could be the third centroid.
 *
 * The order in which the two centroids are passed in /does/ matter wrt spectral-ness. If the first
 * supplied catalog star is `a`, second is `b`, and the unidentified is `c`, then (`a` crossproduct
 * `b`) dot `c` will be positive for all returned `c`. Intuitively, if the triangle is fairly small
 * on the surface of the sphere, then viewed from the origin, the vertices `a`, `b`, and `c` should
 * be clockwise (then, the cross product of `a` and `b` will point mostly towards `c`, and a little
 * outwards, so that dot product with `c` captures that positive outward part).
 */
std::vector<int16_t> IdentifyThirdStar(const PairDistanceKVectorDatabase &db,
                                       const Catalog &catalog,
                                       int16_t catalogIndex1, int16_t catalogIndex2,
                                       float distance1, float distance2,
                                       float tolerance) {

    const int16_t *query1End;
    const int16_t *query1 = db.FindPairsExact(catalog, distance1-tolerance, distance1+tolerance, &query1End);

    const Vec3 &spatial1 = catalog[catalogIndex1].spatial;
    const Vec3 &spatial2 = catalog[catalogIndex2].spatial;
    const Vec3 cross = spatial1.CrossProduct(spatial2);

    // Use PairDistanceInvolvingIterator to find catalog candidates for the unidentified centroid from both sides.

    std::vector<int16_t> result;
    // find all the catalog stars that are in both annuli
    for (PairDistanceInvolvingIterator candidateIt(query1, query1End, catalogIndex1);
         candidateIt.HasValue();
         ++candidateIt) {

        Vec3 candidateSpatial = catalog[*candidateIt].spatial;

        float angle2 = AngleUnit(candidateSpatial, spatial2);

        // check distance to second star
        if (!(angle2 >= distance2-tolerance && angle2 <= distance2+tolerance)) {
            continue;
        }

        // check spectrality
        float spectralTorch = cross * candidateSpatial;
        // if they are nearly coplanar, don't need to check spectrality
        // TODO: Implement ^^. Not high priority, since always checking spectrality is conservative.
        if (spectralTorch <= 0) {
            continue;
        }

        // we've made it through the gauntlet!
        result.push_back(*candidateIt);
    }

    return result;
}

IRUnidentifiedCentroid *SelectNextUnidentifiedCentroid(std::vector<IRUnidentifiedCentroid *> *aboveThresholdCentroids,
                                                      std::vector<IRUnidentifiedCentroid *> *belowThresholdCentroids) {
    if (!belowThresholdCentroids->empty()) {
        auto result = belowThresholdCentroids->back();
        belowThresholdCentroids->pop_back();
        return result;
    }

    // need to find the best in aboveThreshold, if any
    auto bestAboveThreshold = std::min_element(aboveThresholdCentroids->begin(), aboveThresholdCentroids->end(),
        [](const IRUnidentifiedCentroid *a, const IRUnidentifiedCentroid *b) {
            return a->bestAngleFrom90 < b->bestAngleFrom90;
        });

    // 10 is arbitrary; but really it should be less than M_PI_2 when set
    if (bestAboveThreshold != aboveThresholdCentroids->end() && (*bestAboveThreshold)->bestAngleFrom90 < 10) {
        auto result = *bestAboveThreshold;
        aboveThresholdCentroids->erase(bestAboveThreshold);
        return result;
    }

    return NULL;
}

const float kAngleFrom90SoftThreshold = M_PI_4; // TODO: tune this

/**
 * Given some identified stars, attempt to identify the rest.
 *
 * Requires a pair distance database to be present. Iterates through the unidentified centroids in
 * an intelligent order, identifying them one by one.
 */
int IdentifyRemainingStarsPairDistance(StarIdentifiers *identifiers,
                                       const Stars &stars,
                                       const PairDistanceKVectorDatabase &db,
                                       const Catalog &catalog,
                                       const Camera &camera,
                                       float tolerance) {
#ifdef LOST_DEBUG_PERFORMANCE
    auto startTimestamp = std::chrono::steady_clock::now();
#endif
    // initialize all unidentified centroids
    std::vector<IRUnidentifiedCentroid> allUnidentifiedCentroids;
    std::vector<IRUnidentifiedCentroid *> aboveThresholdUnidentifiedCentroids;
    std::vector<IRUnidentifiedCentroid *> belowThresholdUnidentifiedCentroids;
    allUnidentifiedCentroids.reserve(stars.size());
    for (size_t i = 0; i < stars.size(); i++) {
        allUnidentifiedCentroids.push_back(IRUnidentifiedCentroid(stars[i], i));
    }
    // add everything from allUnidentifiedCentroids to above threshold
    aboveThresholdUnidentifiedCentroids.reserve(allUnidentifiedCentroids.size());
    for (size_t i = 0; i < allUnidentifiedCentroids.size(); i++) {
        // only add if index is not equal to any starIndex in identifiers already
        if (std::find_if(identifiers->begin(), identifiers->end(),
            [&allUnidentifiedCentroids, i](const StarIdentifier &identifier) {
                return identifier.starIndex == allUnidentifiedCentroids[i].index;
            }) == identifiers->end()) {

            aboveThresholdUnidentifiedCentroids.push_back(&allUnidentifiedCentroids[i]);
        }
    }

    // sort unidentified centroids by x coordinate
    // Don't need this until we fix the Find thing
    // std::sort(unidentifiedCentroids.begin(), unidentifiedCentroids.end(),
    //     [](const IRUnidentifiedCentroid &a, const IRUnidentifiedCentroid &b) {
    //         return a.star->position.x < b.star->position.x;
    //     });

    // for each identified star, add it to the list of identified stars for each unidentified centroid within range
    for (const auto &starId : *identifiers) {
        AddToAllUnidentifiedCentroids(starId, stars,
                                      &aboveThresholdUnidentifiedCentroids, &belowThresholdUnidentifiedCentroids,
                                      db.MinDistance(), db.MaxDistance(),
                                      kAngleFrom90SoftThreshold,
                                      camera);
    }

    int numExtraIdentifiedStars = 0;

    // keep getting the best unidentified centroid and identifying it
    while (!belowThresholdUnidentifiedCentroids.empty() || !aboveThresholdUnidentifiedCentroids.empty()) {
        IRUnidentifiedCentroid *nextUnidentifiedCentroid
            = SelectNextUnidentifiedCentroid(&aboveThresholdUnidentifiedCentroids, &belowThresholdUnidentifiedCentroids);
        if (nextUnidentifiedCentroid == NULL) {
            break;
        }

        // Project next stars to 3d, find angle between them and current unidentified centroid
        Vec3 unidentifiedSpatial = camera.CameraToSpatial(nextUnidentifiedCentroid->star->position);
        Vec3 spatial1 = camera.CameraToSpatial(stars[nextUnidentifiedCentroid->bestStar1.starIndex].position);
        Vec3 spatial2 = camera.CameraToSpatial(stars[nextUnidentifiedCentroid->bestStar2.starIndex].position);
        float d1 = Angle(spatial1, unidentifiedSpatial);
        float d2 = Angle(spatial2, unidentifiedSpatial);
        float spectralTorch = spatial1.CrossProduct(spatial2) * unidentifiedSpatial;

        // find all the catalog stars that are in both annuli
        // flip arguments for appropriate spectrality.
        std::vector<int16_t> candidates =
            spectralTorch > 0
            ? IdentifyThirdStar(db,
                                catalog,
                                nextUnidentifiedCentroid->bestStar1.catalogIndex,
                                nextUnidentifiedCentroid->bestStar2.catalogIndex,
                                d1, d2, tolerance)
            : IdentifyThirdStar(db,
                                catalog,
                                nextUnidentifiedCentroid->bestStar2.catalogIndex,
                                nextUnidentifiedCentroid->bestStar1.catalogIndex,
                                d2, d1, tolerance);

        if (candidates.size() != 1) { // if there is not exactly one candidate, we can't identify the star. Just remove it from the list.
            if (candidates.size() > 1) {
                std::cerr << "WARNING: Multiple catalog stars matched during identify remaining stars. This should be rare." << std::endl;
            }
        } else {
            // identify the centroid
            identifiers->emplace_back(nextUnidentifiedCentroid->index, candidates[0]);

            // update nearby unidentified centroids with the new identified star
            AddToAllUnidentifiedCentroids(identifiers->back(), stars,
                                          &aboveThresholdUnidentifiedCentroids, &belowThresholdUnidentifiedCentroids,
                                          db.MinDistance(), db.MaxDistance(),
                                          // TODO should probably tune this:
                                          kAngleFrom90SoftThreshold,
                                          camera);

            ++numExtraIdentifiedStars;
        }
    }

    // Select() should always empty out this list
    assert(belowThresholdUnidentifiedCentroids.empty());

#ifdef LOST_DEBUG_PERFORMANCE
    auto endTimestamp = std::chrono::steady_clock::now();
    std::cout << "IdentifyRemainingStarsPairDistance took " << std::chrono::duration_cast<std::chrono::microseconds>(endTimestamp - startTimestamp).count() << "us" << std::endl;
#endif

    return numExtraIdentifiedStars;
}

// allCentroidSpatials should be unit vectors
std::vector<BestPyramidAtStar> ComputeBestPyramids(const std::vector<Vec3> &allCentroidSpatials,
                                                   const std::vector<int16_t> &centroidIndices,
                                                   float minDistance,
                                                   float maxDistance) {
    assert(centroidIndices.size() >= 4);
    assert(allCentroidSpatials.size() >= centroidIndices.size());
    assert(minDistance < maxDistance);
    if (allCentroidSpatials.size() > 0) {
        // the only way this can happen by chance is if it's right on the boresight. Some test will catch it :)
        assert(abs(allCentroidSpatials[0].Magnitude() - 1) < 0.0001);
    }

    std::vector<BestPyramidAtStar> result;
    float maxCos = cos(minDistance);
    float minCos = cos(maxDistance);
    for (int i = 0; i < (int)centroidIndices.size(); i++) {
        // Find the 3 closest centroids to this centroid
        std::vector<std::pair<float, int16_t>> cosines;

        // TODO: optimize this using a sorted centroids list
        for (int j = 0; j < (int)centroidIndices.size(); j++) {
            if (i == j) {
                continue;
            }
            float curCos = allCentroidSpatials[centroidIndices[i]] * allCentroidSpatials[centroidIndices[j]];
            // float curDistance = (allCentroids[centroidIndices[i]].position - allCentroids[centroidIndices[j]].position).Magnitude();
            if (minCos <= curCos && curCos <= maxCos) {
                // emplace the NEGATIVE cosine so that the sort will be smallest angle first.
                cosines.emplace_back(-curCos, centroidIndices[j]);
            }
        }

        if (cosines.size() < 3) {
            // Not enough centroids to make a pyramid
            result.emplace_back(i);
            continue;
        }

        std::partial_sort(cosines.begin(), cosines.begin()+4, cosines.end()); // operator< is defined on pairs to sort by first element

        // Compute the sum of the distances between the 3 closest centroids
        float distancesSum = 0;
        for (int j = 0; j < 3; j++) {
            distancesSum += acos(-cosines[j].first);
        }

        // Add the best pyramid starting from this centroid to the result
        result.emplace_back(centroidIndices[i], cosines[0].second, cosines[1].second, cosines[2].second, distancesSum);
    }

    return result;
}

// see star-id-private.hpp for docs on PyramidIterator
PyramidIterator::PyramidIterator(const std::vector<Vec3> &centroidSpatials, float minDistance, float maxDistance)
    : minDistance(minDistance), maxDistance(maxDistance), allCentroidSpatials(centroidSpatials) {
    for (int i = 0; i < (int)centroidSpatials.size(); i++) {
        untriedCentroidIndices.push_back(i);
    }
}

BestPyramidAtStar PyramidIterator::Next() {
    if (untriedCentroidIndices.size() < 4) {
        return BestPyramidAtStar(-1);
    }

    // Find the best pyramid
    std::vector<BestPyramidAtStar> bestPyramids = ComputeBestPyramids(allCentroidSpatials, untriedCentroidIndices,
                                                                      minDistance, maxDistance);
    assert(!bestPyramids.empty());

    // Find the best pyramid
    auto minIt = std::min_element(bestPyramids.begin(), bestPyramids.end());
    assert(minIt != bestPyramids.end());
    BestPyramidAtStar bestPyramid = *minIt;

    if (bestPyramid.distancesSum < 0) {
        // No suitable pyramid exists
        return BestPyramidAtStar(-1);
    }

    // Remove all the stars in the best pyramid from the list of untried stars
    for (int i = 0; i < 4; i++) {
        // Possible to optimize this using remove_if, or otherwise doing all the removal at once.
        untriedCentroidIndices.erase(std::remove(untriedCentroidIndices.begin(), untriedCentroidIndices.end(),
                                                 bestPyramid.centroidIndices[i]),
                                     untriedCentroidIndices.end());
    }

    return bestPyramid;
}

// IdentifyPyramidResult IdentifyPyramid(const PairDistanceKVectorDatabase &db,
//                                       const Catalog &catalog,
//                                       float tolerance,
//                                       const Vec3 &iSpatial, const Vec3 &jSpatial, const Vec3 &kSpatial, const Vec3 &rSpatial,
//                                       int *a, int *b, int *c, int *d) {

//     // sign of determinant, to detect flipped patterns
//     bool spectralTorch = iSpatial.CrossProduct(jSpatial)*kSpatial > 0;

//     float ijDist = AngleUnit(iSpatial, jSpatial);
//     float ikDist = AngleUnit(iSpatial, kSpatial);
//     float irDist = AngleUnit(iSpatial, rSpatial);
//     float jkDist = AngleUnit(jSpatial, kSpatial);
//     float jrDist = AngleUnit(jSpatial, rSpatial);
//     float krDist = AngleUnit(kSpatial, rSpatial); // TODO: we don't really need to
//     // check krDist, if k has been
//     // verified by i and j it's fine.

//     // we check the distances with the extra tolerance requirement to ensure that
//     // there isn't some pyramid that's just outside the database's bounds, but
//     // within measurement tolerance of the observed pyramid, since that would
//     // possibly cause a non-unique pyramid to be identified as unique.
// #define _CHECK_DISTANCE(_dist) if (_dist < db.MinDistance() + tolerance || _dist > db.MaxDistance() - tolerance) { std::cerr << _dist << std::endl; return IdentifyPyramidResult::CannotPossiblyIdentify; }
//     _CHECK_DISTANCE(ijDist);
//     _CHECK_DISTANCE(ikDist);
//     _CHECK_DISTANCE(irDist);
// #undef _CHECK_DISTANCE

// #if LOST_DEBUG > 2
//     std::cerr << "Passed distance checks for pyramid" << std::endl;
// #endif

//     const int16_t *ijEnd, *ikEnd, *irEnd;
//     const int16_t *const ijQuery = db.FindPairsExact(catalog, ijDist - tolerance, ijDist + tolerance, &ijEnd);
//     const int16_t *const ikQuery = db.FindPairsExact(catalog, ikDist - tolerance, ikDist + tolerance, &ikEnd);
//     const int16_t *const irQuery = db.FindPairsExact(catalog, irDist - tolerance, irDist + tolerance, &irEnd);

//     // TODO: theoretically, it's fastest to relabel j,k,r to put the shortest one first (probably).
//     // But since we're not even sure of that, we'll just leave them alone for now.

// #if LOST_DEBUG > 3
//     std::cerr << "Number of ij pairs: " << (ijEnd-ijQuery)/2 << std::endl;
// #endif
//     std::multimap<int16_t, int16_t> ikMap = PairDistanceQueryToMap(ikQuery, ikEnd);
//     std::multimap<int16_t, int16_t> irMap = PairDistanceQueryToMap(irQuery, irEnd);

//     bool foundMatchYet = false;
//     for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd; iCandidateQuery++) {
//         int iCandidate = *iCandidateQuery;
//         // depending on parity, the first or second star in the pair is the "other" one
//         int jCandidate = (iCandidateQuery - ijQuery) % 2 == 0
//             ? iCandidateQuery[1]
//             : iCandidateQuery[-1];

//         // if debug>4, print the candidates
// #if LOST_DEBUG > 3
//         std::cerr << "iCandidate: " << *iCandidateQuery << std::endl;
//         std::cerr << "jCandidate: " << jCandidate << std::endl;
// #endif

//         const Vec3 &iCandidateSpatial = catalog[iCandidate].spatial;
//         const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;

//         Vec3 ijCandidateCross = iCandidateSpatial.CrossProduct(jCandidateSpatial);

//         for (auto kCandidateIt = ikMap.equal_range(iCandidate); kCandidateIt.first != kCandidateIt.second; kCandidateIt.first++) {
//             // kCandidate.first is iterator, then ->second is the value (other star)
//             int kCandidate = kCandidateIt.first->second;
//             Vec3 kCandidateSpatial = catalog[kCandidate].spatial;
//             bool candidateSpectralTorch = ijCandidateCross*kCandidateSpatial > 0;
//             // checking the spectral-ity early to fail fast
//             if (candidateSpectralTorch != spectralTorch) {
// #if LOST_DEBUG > 3
//                 std::cerr << "skipping candidate " << iCandidate << " " << jCandidate << " " << kCandidate << " because spectral-ity mismatch" << std::endl;
// #endif
//                 continue;
//             }

//             // small optimization: We can calculate jk before iterating through r, so we will!
//             float jkCandidateDist = AngleUnit(jCandidateSpatial, kCandidateSpatial);
//             if (jkCandidateDist < jkDist - tolerance || jkCandidateDist > jkDist + tolerance) {
//                 continue;
//             }

//             // TODO: if there are no jr matches, there's no reason to
//             // continue iterating through all the other k-s. Possibly
//             // enumarete all r matches, according to ir, before this loop
//             // However, on 2023-03-21, I tried this, and according to callgrind under -O3 it didn't improve anything.
//             for (auto rCandidateIt = irMap.equal_range(iCandidate); rCandidateIt.first != rCandidateIt.second; rCandidateIt.first++) {
//                 int rCandidate = rCandidateIt.first->second;
//                 const Vec3 &rCandidateSpatial = catalog[rCandidate].spatial;
//                 float jrCandidateDist = AngleUnit(jCandidateSpatial, rCandidateSpatial);
//                 float krCandidateDist;
//                 if (jrCandidateDist < jrDist - tolerance || jrCandidateDist > jrDist + tolerance) {
//                     continue;
//                 }
//                 krCandidateDist = AngleUnit(kCandidateSpatial, rCandidateSpatial);
//                 if (krCandidateDist < krDist - tolerance || krCandidateDist > krDist + tolerance) {
//                     continue;
//                 }

//                 // we have a match!
//                 if (!foundMatchYet) {
//                     *a = iCandidate;
//                     *b = jCandidate;
//                     *c = kCandidate;
//                     *d = rCandidate;
//                     foundMatchYet = true;
//                 } else {
//                     // not unique!
//                     return IdentifyPyramidResult::MatchedAmbiguously;
//                 }
//             }
//         }
//     }
//     return foundMatchYet
//         ? IdentifyPyramidResult::MatchedUniquely
//         : IdentifyPyramidResult::NoMatch;

// }

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL || stars.size() < 4) {
        std::cerr << "Not enough stars, or database missing." << std::endl;
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);
    std::vector<Vec3> starSpatials;
    for (const Star &star : stars) {
        starSpatials.push_back(camera.CameraToSpatial(star.position).Normalize());
    }

    // smallest normal single-precision float is around 10^-38 so we should be all good. See
    // Analytic_Star_Pattern_Probability on the HSL wiki for details.
    float expectedMismatchesConstant = pow(numFalseStars, 4) * pow(tolerance, 5) / 2 / pow(M_PI, 2);

    // keep iterating through pyramids
    PyramidIterator pyramidIterator(starSpatials,
                                    vectorDatabase.MinDistance() + tolerance,
                                    vectorDatabase.MaxDistance() - tolerance);
    BestPyramidAtStar bestPyramid(-1);
    while (bestPyramid = pyramidIterator.Next(), !bestPyramid.isNull()) {

#if LOST_DEBUG > 2
        // print distance sum
        std::cerr << "Current pyramid distance sum: " << bestPyramid.distancesSum << std::endl;
        // print out each centroid index
        std::cerr << bestPyramid.centroidIndices[0] << " " << bestPyramid.centroidIndices[1] << " " << bestPyramid.centroidIndices[2] << " " << bestPyramid.centroidIndices[3] << std::endl;
#endif
        int i = bestPyramid.centroidIndices[0],
            j = bestPyramid.centroidIndices[1],
            k = bestPyramid.centroidIndices[2],
            r = bestPyramid.centroidIndices[3];

        assert(i != j && j != k && k != r && i != k && i != r && j != r);

        // TODO: move this out of the loop?
        Vec3 iSpatial = starSpatials[i];
        Vec3 jSpatial = starSpatials[j];
        Vec3 kSpatial = starSpatials[k];

        float ijDist = AngleUnit(iSpatial, jSpatial);

        float iSinInner = sin(Angle(jSpatial - iSpatial, kSpatial - iSpatial));
        float jSinInner = sin(Angle(iSpatial - jSpatial, kSpatial - jSpatial));
        float kSinInner = sin(Angle(iSpatial - kSpatial, jSpatial - kSpatial));

        // if we made it this far, all 6 angles are confirmed! Now check
        // that this match would not often occur due to chance.
        // See Analytic_Star_Pattern_Probability on the HSL wiki for details
        float expectedMismatches = expectedMismatchesConstant
            * sin(ijDist)
            / kSinInner
            / std::max(std::max(iSinInner, jSinInner), kSinInner);

        if (expectedMismatches > maxMismatchProbability) {
            std::cout << "skip: mismatch prob." << std::endl;
            continue;
        }

        Vec3 rSpatial = starSpatials[r];

        int iMatch, jMatch, kMatch, rMatch;
        // int numMatches = identifyResult
        //     = IdentifyPyramid(vectorDatabase, catalog, tolerance,
        //                       iSpatial, jSpatial, kSpatial, rSpatial,
        //                       &iMatch, &jMatch, &kMatch, &rMatch);

        // switch (identifyResult) {
        // case IdentifyPyramidResult::CannotPossiblyIdentify:
        // case IdentifyPyramidResult::NoMatch:
        // case IdentifyPyramidResult::MatchedAmbiguously:
        //     continue;
        // case IdentifyPyramidResult::MatchedUniquely:
        //     identified.push_back(StarIdentifier(i, iMatch));
        //     identified.push_back(StarIdentifier(j, jMatch));
        //     identified.push_back(StarIdentifier(k, kMatch));
        //     identified.push_back(StarIdentifier(r, rMatch));

        //     int numAdditionallyIdentified = IdentifyRemainingStarsPairDistance(&identified, stars, vectorDatabase, catalog, camera, tolerance);
        //     printf("Identified an additional %d stars.\n", numAdditionallyIdentified);
        //     assert(numAdditionallyIdentified == (int)identified.size()-4);

        //     return identified;
        // }
    }

    std::cerr << "Tried all pyramids; none matched." << std::endl;
    return identified;
}

}
