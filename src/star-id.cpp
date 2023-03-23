#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <chrono>
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
    std::multimap<int16_t, int16_t> result;
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

    // smallest normal single-precision float is around 10^-38 so we should be all good. See
    // Analytic_Star_Pattern_Probability on the HSL wiki for details.
    float expectedMismatchesConstant = pow(numFalseStars, 4) * pow(tolerance, 5) / 2 / pow(M_PI, 2);

    // this iteration technique is described in the Pyramid paper. Briefly: i will always be the
    // lowest index, then dj and dk are how many indexes ahead the j-th star is from the i-th, and
    // k-th from the j-th. In addition, we here add some other numbers so that the pyramids are not
    // weird lines in wide FOV images. TODO: Select the starting points to ensure that the first pyramids are all within measurement tolerance.
    int numStars = (int)stars.size();
    // the idea is that the square root is about across the FOV horizontally
    int across = floor(sqrt(numStars))*2;
    int halfwayAcross = floor(sqrt(numStars)/2);
    long totalIterations = 0;

    int jMax = numStars - 3;
    for (int jIter = 0; jIter < jMax; jIter++) {
        int dj = 1+(jIter+halfwayAcross)%jMax;

        int kMax = numStars-dj-2;
        for (int kIter = 0; kIter < kMax; kIter++) {
            int dk = 1+(kIter+across)%kMax;

            int rMax = numStars-dj-dk-1;
            for (int rIter = 0; rIter < rMax; rIter++) {
                int dr = 1+(rIter+halfwayAcross)%rMax;

                int iMax = numStars-dj-dk-dr-1;
                for (int iIter = 0; iIter <= iMax; iIter++) {
                    int i = (iIter + iMax/2)%(iMax+1); // start near the center of the photo

                    // identification failure due to cutoff
                    if (++totalIterations > cutoff) {
                        std::cerr << "Cutoff reached." << std::endl;
                        return identified;
                    }

                    int j = i+dj;
                    int k = j+dk;
                    int r = k+dr;

                    assert(i != j && j != k && k != r && i != k && i != r && j != r);

                    // TODO: move this out of the loop?
                    Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
                    Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
                    Vec3 kSpatial = camera.CameraToSpatial(stars[k].position).Normalize();

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

                    Vec3 rSpatial = camera.CameraToSpatial(stars[r].position).Normalize();

                    // sign of determinant, to detect flipped patterns
                    bool spectralTorch = iSpatial.CrossProduct(jSpatial)*kSpatial > 0;

                    float ikDist = AngleUnit(iSpatial, kSpatial);
                    float irDist = AngleUnit(iSpatial, rSpatial);
                    float jkDist = AngleUnit(jSpatial, kSpatial);
                    float jrDist = AngleUnit(jSpatial, rSpatial);
                    float krDist = AngleUnit(kSpatial, rSpatial); // TODO: we don't really need to
                                                                  // check krDist, if k has been
                                                                  // verified by i and j it's fine.

                    // we check the distances with the extra tolerance requirement to ensure that
                    // there isn't some pyramid that's just outside the database's bounds, but
                    // within measurement tolerance of the observed pyramid, since that would
                    // possibly cause a non-unique pyramid to be identified as unique.
#define _CHECK_DISTANCE(_dist) if (_dist < vectorDatabase.MinDistance() + tolerance || _dist > vectorDatabase.MaxDistance() - tolerance) { continue; }
                    _CHECK_DISTANCE(ikDist);
                    _CHECK_DISTANCE(irDist);
                    _CHECK_DISTANCE(jkDist);
                    _CHECK_DISTANCE(jrDist);
                    _CHECK_DISTANCE(krDist);
#undef _CHECK_DISTANCE

                    const int16_t *ijEnd, *ikEnd, *irEnd;
                    const int16_t *const ijQuery = vectorDatabase.FindPairsLiberal(ijDist - tolerance, ijDist + tolerance, &ijEnd);
                    const int16_t *const ikQuery = vectorDatabase.FindPairsLiberal(ikDist - tolerance, ikDist + tolerance, &ikEnd);
                    const int16_t *const irQuery = vectorDatabase.FindPairsLiberal(irDist - tolerance, irDist + tolerance, &irEnd);

                    std::multimap<int16_t, int16_t> ikMap = PairDistanceQueryToMap(ikQuery, ikEnd);
                    std::multimap<int16_t, int16_t> irMap = PairDistanceQueryToMap(irQuery, irEnd);

                    int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
                    for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd; iCandidateQuery++) {
                        int iCandidate = *iCandidateQuery;
                        // depending on parity, the first or second star in the pair is the "other" one
                        int jCandidate = (iCandidateQuery - ijQuery) % 2 == 0
                            ? iCandidateQuery[1]
                            : iCandidateQuery[-1];

                        const Vec3 &iCandidateSpatial = catalog[iCandidate].spatial;
                        const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;

                        Vec3 ijCandidateCross = iCandidateSpatial.CrossProduct(jCandidateSpatial);

                        for (auto kCandidateIt = ikMap.equal_range(iCandidate); kCandidateIt.first != kCandidateIt.second; kCandidateIt.first++) {
                            // kCandidate.first is iterator, then ->second is the value (other star)
                            int kCandidate = kCandidateIt.first->second;
                            Vec3 kCandidateSpatial = catalog[kCandidate].spatial;
                            bool candidateSpectralTorch = ijCandidateCross*kCandidateSpatial > 0;
                            // checking the spectral-ity early to fail fast
                            if (candidateSpectralTorch != spectralTorch) {
                                continue;
                            }

                            // small optimization: We can calculate jk before iterating through r, so we will!
                            float jkCandidateDist = AngleUnit(jCandidateSpatial, kCandidateSpatial);
                            if (jkCandidateDist < jkDist - tolerance || jkCandidateDist > jkDist + tolerance) {
                                continue;
                            }

                            // TODO: if there are no jr matches, there's no reason to
                            // continue iterating through all the other k-s. Possibly
                            // enumarete all r matches, according to ir, before this loop
                            for (auto rCandidateIt = irMap.equal_range(iCandidate); rCandidateIt.first != rCandidateIt.second; rCandidateIt.first++) {
                                int rCandidate = rCandidateIt.first->second;
                                const Vec3 &rCandidateSpatial = catalog[rCandidate].spatial;
                                float jrCandidateDist = AngleUnit(jCandidateSpatial, rCandidateSpatial);
                                float krCandidateDist;
                                if (jrCandidateDist < jrDist - tolerance || jrCandidateDist > jrDist + tolerance) {
                                    continue;
                                }
                                krCandidateDist = AngleUnit(kCandidateSpatial, rCandidateSpatial);
                                if (krCandidateDist < krDist - tolerance || krCandidateDist > krDist + tolerance) {
                                    continue;
                                }

                                // we have a match!

                                if (iMatch == -1) {
                                    iMatch = iCandidate;
                                    jMatch = jCandidate;
                                    kMatch = kCandidate;
                                    rMatch = rCandidate;
                                } else {
                                    // uh-oh, stinky!
                                    // TODO: test duplicate detection, it's hard to cause it in the real catalog...
                                    std::cerr << "Pyramid not unique, skipping..." << std::endl;
                                    goto sensorContinue;
                                }
                            }
                        }


                    }

                    if (iMatch != -1) {
                        printf("Matched unique pyramid!\nExpected mismatches: %e\n", expectedMismatches);
                        identified.push_back(StarIdentifier(i, iMatch));
                        identified.push_back(StarIdentifier(j, jMatch));
                        identified.push_back(StarIdentifier(k, kMatch));
                        identified.push_back(StarIdentifier(r, rMatch));

                        int numAdditionallyIdentified = IdentifyRemainingStarsPairDistance(&identified, stars, vectorDatabase, catalog, camera, tolerance);
                        printf("Identified an additional %d stars.\n", numAdditionallyIdentified);
                        assert(numAdditionallyIdentified == (int)identified.size()-4);

                        return identified;
                    }

                sensorContinue:;
                }
            }
        }
    }

    std::cerr << "Tried all pyramids; none matched." << std::endl;
    return identified;
}

/**
 * Try to identify a pattern of starts using pair distances.
 *
 * Pass it 3D spatials, and it will try to match catalog indices
 * @return the "number" of matches. -1 if it cannot possibly match because the inter-star distances are outside the allowable range. 0 if just no match. 1 if unique match (this is the only case you can rely on the values of matchedCatalogIndices). Or, if not a unique match, the number of matches.
 */
template <int numPatternStars>
std::vector<std::vector<int16_t>> IdentifyPatternPairDistance(const PairDistanceKVectorDatabase &db,
                                                                  const Catalog &catalog,
                                                                  float tolerance,
                                                                  const Vec3 spatials[numPatternStars]) {
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
                return {};
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

    std::vector<std::vector<int16_t>> result;

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

            {
                std::vector<int16_t> currentMatch;
                std::copy(candidateCatalogIndices, candidateCatalogIndices + numPatternStars, std::back_inserter(currentMatch));
                result.push_back(currentMatch);
            }
            // now draw the rest of the fucking owl. Idea: For every remaining star, find stars that
            // distance away from the 0-th star using a hashmap (just like we did for the 2nd star),
            // then verify its distances to all other stars.
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

                    // check distances against all other stars, except 0, because that's part of the query hence already checked
                    for (int m = 1; m < l; m++) {
                        float lmCandidateDist = AngleUnit(candidateSpatials[l], candidateSpatials[m]);
                        if (lmCandidateDist < distances[l][m] - tolerance || lmCandidateDist > distances[l][m] + tolerance) {
                            goto nextLCandidate;
                        }
                    }

                    // if there are multiple matches, the matched indices are undefined anyway, so just set them unconditionally

                nextLCandidate:;
                }
            }
        }
    nextICandidate:;
    }

    return result;
}

BayesianStarIdAlgorithm::BayesianStarIdAlgorithm(
    float tolerance, int numFalseStars,
    float softConfidenceThreshold, float hardConfidenceThreshold,
    float admissibleIgnoredProbability):
    tolerance(tolerance), numFalseStars(numFalseStars),
    softConfidenceThreshold(softConfidenceThreshold), hardConfidenceThreshold(hardConfidenceThreshold),
    admissibleIgnoredProbability(admissibleIgnoredProbability) {

    assert(0.5 <= hardConfidenceThreshold); // some of our algorithm uses this
    assert(hardConfidenceThreshold <= softConfidenceThreshold);
    assert(softConfidenceThreshold <= 1.0);
}

class BayesPossibility {
public:
    BayesPossibility(float probability)
        : probability(probability) { };
    BayesPossibility(float probability, int16_t centroidIndex)
        : probability(probability) {
        centroidIndices.push_back(centroidIndex);
    };
    BayesPossibility(float probability, int16_t centroidIndex1, int16_t centroidIndex2,
                     int numPairs)
        : probability(probability), numPairs(numPairs) {
        centroidIndices.push_back(centroidIndex1);
        centroidIndices.push_back(centroidIndex2);
    }
    BayesPossibility() = default;

    // the number of subconfigurations stored in this possibility (subconfigurations are really just
    // more possibilities)
    int NumConfigurations(const Catalog &catalog) const {
        switch (NumTrueStars()) {
        case 2:
            return numPairs;
        case 1:
            return catalog.size();
        case 0:
            return 1;
        default:
            assert(!catalogIndices.empty());
            return catalogIndices.size() / NumTrueStars();
        }
    }

    int NumTrueStars() const {
        return centroidIndices.size();
    }

    float TotalProbability(const Catalog &catalog) const {
        return probability * NumConfigurations(catalog);
    }

    // convert this possibility to an identifier if appropriate
    StarIdentifiers ToStarIdentifiers() const {
        assert(centroidIndices.size() > 2);
        // no NumConfigurations() just so we avoid needing catalog
        assert(centroidIndices.size() == catalogIndices.size());
        StarIdentifiers result;
        for (int i = 0; i < (int)centroidIndices.size(); i++) {
            result.push_back(StarIdentifier(centroidIndices[i], catalogIndices[i]));
        }
        return result;
    };

    float probability; // Probability of each possible configuration inside
    std::vector<int16_t> centroidIndices;
    int numPairs; // used when NumTrueStars()==2
    std::vector<int16_t> catalogIndices;
};

// // can iterate over either manually passed in int16_t*s, or vector iterators. When iterating over
// // int16_t*s, also goes over them in reverse
// class BayesPossibilityIterator {
// public:
//     BayesPossibilityIterator(BayesPossibility possibility)
//         : queryStart(possibility.catalogIndicesBegin), queryCur(possibility.catalogIndicesBegin), queryEnd(possibility.catalogIndicesEnd),
//           vecCur(possibility.catalogIndices.cbegin()), vecEnd(possibility.catalogIndices.cend())
//         { };

//     bool hasValue() const {
//         if (queryCur != NULL) {
//             return !(queryCur == queryStart && reverse == true);
//         } else {
//             return vecCur != vecEnd;
//         }
//     }
//     void operator++() {
//         assert(hasValue());
//         if (queryCur != NULL) {
//             if (!reverse) {
//                 queryCur++;
//                 if (queryCur == queryEnd) {
//                     queryCur--;
//                     reverse = true;
//                 }
//             } else {
//                 queryCur--;
//             }
//         } else {
//             vecCur++;
//         }
//     };
//     int16_t operator*() {
//         assert(hasValue());
//         if (queryCur != NULL) {
//             return *queryCur;
//         } else {
//             return *vecCur;
//         }
//     };

// private:
//     const int16_t *queryStart;
//     const int16_t *queryCur = NULL;
//     const int16_t *queryEnd;
//     bool reverse = false;
//     std::vector<int16_t>::const_iterator vecCur;
//     std::vector<int16_t>::const_iterator vecEnd;
// };

typedef std::vector<BayesPossibility> BayesPrior;

// summed over all possibilities
static float BayesPriorTotalProbability(const BayesPrior &prior, const Catalog &catalog) {
    float result = 0;
    for (const BayesPossibility &possibility : prior) {
        result += possibility.TotalProbability(catalog);
    }
    return result;
}

class DebugBayesPossibilitiesSummary;
typedef std::vector<DebugBayesPossibilitiesSummary> DebugBayesPriorSummary;
static DebugBayesPriorSummary DebugCalculateBayesPriorSummary(const BayesPrior &prior, const Catalog &catalog);
static void DebugPrintBayesPriorSummary(const DebugBayesPriorSummary &summary);

// sorts possibilities by /per-configuration/ probability decreasing
static bool PossibilityProbabilityComparator(const BayesPossibility &p1, const BayesPossibility &p2) {
    return p1.probability > p2.probability;
}

// return a set of star identfiers that has confidence above the given threshold for all stars. Will
// sort the prior in-place
static StarIdentifiers Mode(BayesPrior *prior, const Catalog &catalog, float hardConfidenceThreshold, int numCentroids) {
    std::sort(prior->begin(), prior->end(), &PossibilityProbabilityComparator);
    // this could be passed in as an argument, since it's already calculated, but putting it here
    // avoids bugs.
    float totalProbability = BayesPriorTotalProbability(*prior, catalog);

    // TODO: a lot of stuff in bayes should be moved to int16_t to save memory
    std::vector<int16_t> centroidIds(numCentroids, -1);
    std::vector<float> centroidProbs(numCentroids, 0.0);
    for (const BayesPossibility &possibility : *prior) {
        // for now, assume that only possibilities with >= 3 stars are going to help us.
        // TODO: think about this assumption. Can 2-star possibilities really help us?
        if (possibility.NumTrueStars() < 3) {
            continue;
        }

        for (auto catalogIndexIt = possibility.catalogIndices.begin();
             catalogIndexIt != possibility.catalogIndices.end();) {

            for (int i = 0; i < possibility.NumTrueStars(); i++, catalogIndexIt++) {
                int centroidIndex = possibility.centroidIndices[i];
                if (centroidIds[centroidIndex] == -1
                    && possibility.probability/totalProbability >= (1-hardConfidenceThreshold)) {

                    centroidIds[centroidIndex] = *catalogIndexIt;
                }

                if (*catalogIndexIt == centroidIds[centroidIndex]) {
                    centroidProbs[centroidIndex] += possibility.probability/totalProbability;
                }
            }
        }
    }

    StarIdentifiers result;
    for (int i = 0; i < (int)centroidIds.size(); i++) {
        if (centroidIds[i] != -1 && centroidProbs[i] >= hardConfidenceThreshold) {
            result.push_back(StarIdentifier(i, centroidIds[i]));
        }
    }
    return result;

    // float sum = 0.0;
    // float max = -1.0;
    // BayesPossibility maxPossibility;
    // for (const BayesPossibility &possibility : prior) {
    //     sum += possibility.TotalProbability(catalog);
    //     // TODO: check that numconfigurations always storable in 16 bits? Hopefully should be, but
    //     // maybe just assert?
    //     int numConfigurations = possibility.NumConfigurations(catalog);
    //     if (numConfigurations == 1 && possibility.probability > max) {
    //         max = possibility.probability;
    //         maxPossibility = possibility;
    //     }
    // }
    // if (max < 0 || maxPossibility.NumTrueStars() <= 2) {
    //     *modeProbability = -1.0;
    //     return StarIdentifiers();
    // } else {
    //     *modeProbability = max/sum;
    //     return maxPossibility.ToStarIdentifiers();
    // }
}

// find approximate area of intersection between 
static float AnnulusIntersectionArea(float tolerance,
                                     const Vec3 &knownStar1, const Vec3 &knownStar2,
                                     const Vec3 &newStar) {
    const Vec3 diff1 = (knownStar1-newStar).Normalize();
    const Vec3 diff2 = (knownStar2-newStar).Normalize();
    // could also just do sin(AngleUnit(diff1, diff2)), but it's fancier to use trig identities
    const float sinInnerAngle = sqrt(1 - diff1*diff2);
    return 4*tolerance*tolerance/sinInnerAngle;
}

StarIdentifiers BayesianStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars,
                                            const Catalog &catalog, const Camera &camera) const {
    const float kUnitSphereArea = 4*M_PI;
    // for each new centroid, we consider information from up to kMaxNearbyCentroids many centroids
    const int kIdealNumNeighbors = 4;
    const int kMinNumNeighbors = 2;
    // // if there are fewer than this many possibilities for a star, check them all individually
    // // rather than using an involving iterator to narrow them down.
    // const int kInvolvingIteratorCutoff = 30;

    // load database
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        std::cerr << "Database missing" << std::endl;
        return StarIdentifiers();
    }
    PairDistanceKVectorDatabase kvector(databaseBuffer);

    // normalized spatials for each centroid
    std::vector<Vec3> starSpatials;
    for (const Star &star : stars) {
        starSpatials.push_back(camera.CameraToSpatial(star.position).Normalize());
    }

    // we only keep track of un-normalized probabilities (do not sum to one). So, we can multiply
    // all of them by an arbitrary constant without functionally changing anything. This variable is
    // picked just to keep probabilities in a reasonable range. Towards this end, we set constant
    // approximately so that multiplying by 
    float probScale = tolerance*tolerance*3*M_PI;

    // // Probability of a given star occurring at a given point in space, within tolerance. sigma^2/4
    // // comes from area of a circle with radius sigma divided by surface area of a sphere.
    // float specificStarHereLikelihood = tolerance*tolerance/4;
    // // strictly speaking, this is the expected number of false stars in a given location. But this
    // // expected value should be very close to 0, so it's about equal to probability.
    // float falseStarHereLikelihood = numFalseStars * specificStarHereLikelihood;

    // initialize the prior
    BayesPrior prior;
    BayesPossibility noTrueStars(1.0e12);
    prior.push_back(noTrueStars);

    const int firstCentroid = stars.size() / 2; // always strictly <stars.size()
    // we want to consider all centroids, but in an order such that each explored centroid is close
    // to a maximal number of already considered centroids, up to kIdealNumNeighbors, beyond which
    // we don't care. To achieve this, repeatedly loop through all centroids. If we loop through and
    // don't find anything with 4 neighbors, switch to "dire straits" mode where we accept with only
    // two neighbors. If still nothing, bail out. An improvement would handle the case when there
    // are multiple groups of centroids which aren't close enough to interact.
    std::vector<int> exploredCentroids;
    enum class CentroidSearchMode {
        NoneFound, // Require min(4, numCentroidsExplored) neighbors, no centroids explored
                   // this round yet.
        Found, // Require min(4, numCentroidsExplored) neighbors, at least one centroid explored
        DireStraits, // Require min(2, numCentroidsExplored) neighbors, no centroids explored 
    };
    CentroidSearchMode searchMode = CentroidSearchMode::Found;
    while (searchMode != CentroidSearchMode::DireStraits) {
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
        const int numCentroidsExplored = exploredCentroids.size();
        int minNumNeighbors;
        if (searchMode == CentroidSearchMode::Found) {
            minNumNeighbors = std::min(kIdealNumNeighbors, numCentroidsExplored);
            // first few iterations are automatically dire straits
            if (minNumNeighbors < kMinNumNeighbors) {
                searchMode = CentroidSearchMode::DireStraits;
            } else {
                searchMode = CentroidSearchMode::NoneFound;
            }
        } else {
            assert(searchMode == CentroidSearchMode::NoneFound);
            minNumNeighbors = std::min(kMinNumNeighbors, numCentroidsExplored);
            searchMode = CentroidSearchMode::DireStraits;
        }

        // TODO: to make it faster to look through centroids, calculate the min/max x/y values to
        // consider, then we can fail fast.

        int curCentroid, numNeighbors;
        bool foundCentroid = false;
        for (int centroidOffset = 0; centroidOffset < (int)stars.size(); centroidOffset++) {
            curCentroid = (firstCentroid + centroidOffset) % stars.size();

            // ensure that there are enough neighbors
            numNeighbors = 0;
            bool alreadyExplored = false;
            // performance improvement could be implemented: Check first that the current centroid
            // is not already explored before even starting this loop.
            for (int exploredCentroid : exploredCentroids) {
                if (exploredCentroid == curCentroid) {
                    alreadyExplored = true;
                    break;
                }

                float distance = AngleUnit(starSpatials[exploredCentroid], starSpatials[curCentroid]);
                if (kvector.DistanceInRange(distance, tolerance)) {
                    numNeighbors++;
                }
            }
            if (alreadyExplored || numNeighbors < minNumNeighbors) {
                continue;
            }
            // There are enough neighbors!
            foundCentroid = true;
            break;
        }

        if (!foundCentroid) {
            continue;
        }

        searchMode = CentroidSearchMode::Found;

        int origPriorLength = prior.size();
        for (int priorI = 0; priorI < origPriorLength; priorI++) {
            // safer to copy than make a reference, because prior.push_back can reallocate it
            BayesPossibility possibility = prior[priorI];
            // // a more theoretically accurate way would not include possibilities where the current
            // // star must be true. Ie, any configuration that indicates a true star at the current
            // // location certainly shouldn't be included. TODO wouldn't be hard to improve
            // priorSum += possibility.TotalProbability(catalog);

            switch (possibility.NumTrueStars()) {
            case 0: {
                float trueProb = possibility.probability/kUnitSphereArea * probScale;
                prior.push_back(BayesPossibility(trueProb, curCentroid));

                break;
            }

            case 1: {
                float distance = AngleUnit(starSpatials[possibility.centroidIndices[0]],
                                           starSpatials[curCentroid]);
                if (!kvector.DistanceInRange(distance, tolerance)) {
                    std::cerr << "WARNING: centroid out of range of a one-true-star possibility!" << std::endl;
                    break;
                    // presently, still insert the possibility with the new star false. Should
                    // look into if we can do something better easily.
                }
                const int16_t *end;
                const int16_t *begin = kvector.FindPairsLiberal(distance-tolerance, distance+tolerance, &end);
                // The number of pairs, in the order returned, is (end-begin)/2. But we want the
                // number of pairs considering both possible orientations, so just end-begin
                int numPairs = end-begin;
                if (begin != end) {
                    // area of spherical annulus with radius sigma
                    const float annulusArea = 2*M_PI*tolerance*sin(distance);
                    const float trueProb = possibility.probability/annulusArea * probScale;
                    prior.push_back(
                        BayesPossibility(trueProb, possibility.centroidIndices[0], curCentroid, numPairs));
                }
                break;
            }

            default: { // >= 2
                ///// STEP 1: Find the closest centroidIndices in the possibility /////
                // unfortunately, this means storing indexes into an indexes array...
                std::vector<int> closestCentroidIndicesIndices;
                // can also do this with std::iota
                for (int i = 0; i < (int)possibility.NumTrueStars(); i++) {
                    closestCentroidIndicesIndices.push_back(i);
                }

                int possibilityNumNeighbors = std::min(minNumNeighbors, possibility.NumTrueStars());

                // find the closest numNeighbors neighbors. Closest will make the rest of this
                // faster and tend to be well distributed around different sides of the new
                // centroid, which keeps the annulus intersection area small.
                // nth_element partitions the array up to the given mid iterator. Still need to sort afterward
                // brrt brrt style guideline violation: no lambdas
                std::nth_element(closestCentroidIndicesIndices.begin(), closestCentroidIndicesIndices.begin() + possibilityNumNeighbors - 1, closestCentroidIndicesIndices.end(),
                                 [curCentroid, &stars, &possibility](int s1, int s2) -> bool {
                                     const Vec2 &p = stars[curCentroid].position;
                                     const Vec2 &p1 = stars[possibility.centroidIndices[s1]].position;
                                     const Vec2 &p2 = stars[possibility.centroidIndices[s2]].position;
                                     const float s1Dist = (p-p1).MagnitudeSq();
                                     const float s2Dist = (p-p2).MagnitudeSq();
                                     return s1Dist < s2Dist;
                                 });
                // TODO: deduplicate
                std::sort(closestCentroidIndicesIndices.begin(), closestCentroidIndicesIndices.begin() + possibilityNumNeighbors,
                                 [curCentroid, &stars, &possibility](int s1, int s2) -> bool {
                                     const Vec2 &p = stars[curCentroid].position;
                                     const Vec2 &p1 = stars[possibility.centroidIndices[s1]].position;
                                     const Vec2 &p2 = stars[possibility.centroidIndices[s2]].position;
                                     const float s1Dist = (p-p1).MagnitudeSq();
                                     const float s2Dist = (p-p2).MagnitudeSq();
                                     return s1Dist < s2Dist;
                                 });

                auto firstNeighborIndexIndex = closestCentroidIndicesIndices.end();
                auto lastNeighborIndexIndex = closestCentroidIndicesIndices.end();
                for (auto centroidIndexIt = closestCentroidIndicesIndices.begin(); centroidIndexIt != closestCentroidIndicesIndices.end(); centroidIndexIt++) {
                    float dist = AngleUnit(starSpatials[curCentroid], starSpatials[possibility.centroidIndices[*centroidIndexIt]]);

                    if (dist > kvector.MaxDistance()) {
                        break;
                    }

                    if (dist >= kvector.MinDistance()) {
                        if (firstNeighborIndexIndex == closestCentroidIndicesIndices.end()) {
                            firstNeighborIndexIndex = centroidIndexIt;
                        }
                        if (dist <= kvector.MaxDistance()) {
                            lastNeighborIndexIndex = centroidIndexIt;
                        }
                    }
                }
                // number of neighbors in this possibility specifically
                int numNeighborsPossibility = lastNeighborIndexIndex-firstNeighborIndexIndex+1;

                if (numNeighborsPossibility < 2) {
                    // TODO: is this a big deal?
                    std::cerr << "WARNING: Not enough neighboring stars in current possibility. It will fall to the wayside. Try increasing database max distance." << std::endl;
                    break;
                }

                BayesPossibility newPossibility;
                newPossibility.centroidIndices = possibility.centroidIndices;
                newPossibility.centroidIndices.push_back(curCentroid);

                // for each pair of stars we're going to consider, estimate the area where the
                // new star is constrained to as the smallest intersection of annuluses around
                // two of the neighbors.
                float smallestArea = INFINITY;
                for (auto neighbor1It = firstNeighborIndexIndex; neighbor1It <= lastNeighborIndexIndex-1; neighbor1It++) {
                    for (auto neighbor2It = firstNeighborIndexIndex+1; neighbor2It <= lastNeighborIndexIndex; neighbor2It++) {
                        float curArea = AnnulusIntersectionArea(
                            tolerance,
                            starSpatials[possibility.centroidIndices[*neighbor1It]],
                            starSpatials[possibility.centroidIndices[*neighbor2It]],
                            starSpatials[curCentroid]);
                        smallestArea = std::min(smallestArea, curArea);
                    }
                }
                newPossibility.probability = possibility.probability/smallestArea * probScale;

                if (possibility.NumTrueStars() == 2) {
                    const Vec3 spatials[3] = {
                        starSpatials[possibility.centroidIndices[*firstNeighborIndexIndex]],
                        starSpatials[possibility.centroidIndices[*lastNeighborIndexIndex]],
                        starSpatials[curCentroid]
                    };
                    auto matchedConfigurations = IdentifyPatternPairDistance<3>(kvector, catalog, tolerance, spatials);
                    std::cerr << matchedConfigurations.size() << std::endl;
                    for (auto &matchedConfiguration : matchedConfigurations) {
                        for (int16_t catalogIndex : matchedConfiguration) {
                            newPossibility.catalogIndices.push_back(catalogIndex);
                        }
                    }
                    // // i is the new star, j and k are the existing stars.
                    // const float ijDist = AngleUnit(starSpatials[curCentroid],
                    //                                starSpatials[possibility.centroidIndices[*firstNeighborIndexIndex]]);
                    // const float ikDist = AngleUnit(starSpatials[curCentroid],
                    //                                starSpatials[possibility.centroidIndices[*lastNeighborIndexIndex]]);
                    // const float jkDist = AngleUnit(starSpatials[possibility.centroidIndices[*firstNeighborIndexIndex]],
                    //                                starSpatials[possibility.centroidIndices[*lastNeighborIndexIndex]]);
                    // const int16_t *ijEnd, *ikEnd;
                    // const int16_t *const ijQuery = kvector.FindPairsLiberal(ijDist - tolerance, ijDist + tolerance, &ijEnd);
                    // const int16_t *const ikQuery = kvector.FindPairsLiberal(ikDist - tolerance, ikDist + tolerance, &ikEnd);

                        
                    // std::vector<bool> iSeen(catalog.size(), false);
                    // for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd; iCandidateQuery++) {
                    //     int iCandidate = *iCandidateQuery;
                    //     if (iSeen[iCandidate]) {
                    //         continue;
                    //     }
                    //     iSeen[iCandidate] = true;

                    //     PairDistanceInvolvingIterator jIterator(ijQuery, ijEnd, iCandidate);
                    //     PairDistanceInvolvingIterator kIterator(ikQuery, ikEnd, iCandidate);
                    //     std::vector<int16_t> jCandidates;
                    //     std::vector<int16_t> kCandidates;
                    //     while (jIterator.HasValue()) {
                    //         jCandidates.push_back(*jIterator);
                    //         ++jIterator;
                    //     }
                    //     while (kIterator.HasValue()) {
                    //         kCandidates.push_back(*kIterator);
                    //         ++kIterator;
                    //     }
                    //     for (int16_t jCandidate : jCandidates) {
                    //         for (int16_t kCandidate : kCandidates) {
                    //             const float jkCatalogDist = abs(AngleUnit(catalog[jCandidate].spatial, catalog[kCandidate].spatial));
                    //             if (abs(jkCatalogDist - jkDist) <= tolerance) {
                    //                 // TODO: check spectral

                    //                 // this is a match! Add the possibility
                    //                 newPossibility.catalogIndices.push_back(jCandidate);
                    //                 newPossibility.catalogIndices.push_back(kCandidate);
                    //                 newPossibility.catalogIndices.push_back(iCandidate);
                    //             }
                    //         }
                    //     }
                    // }
                } else {
                    // >2
                    assert(possibility.NumTrueStars() > 2);

                    const float ijDist = AngleUnit(starSpatials[curCentroid],
                                                   starSpatials[possibility.centroidIndices[*firstNeighborIndexIndex]]);
                    const int16_t *ijEnd;
                    const int16_t *const ijQuery = kvector.FindPairsLiberal(ijDist-tolerance, ijDist+tolerance, &ijEnd);
                    // TODO: it still may be faster to use an involving iterator if there are only a small number of configurations, just because building the map is expensive.
                    // std::multimap<int16_t, int16_t> ijMap = PairDistanceQueryToMap(ijQuery, ijEnd);
                    std::vector<int> oldConfigurationsToRemove;
                    for (int conf = 0; conf < possibility.NumConfigurations(catalog); conf++) {
                        // first, narrow down to stars that are compatible with the first star
                        // in the configuration. Should only be a handful.
                        auto catalogIndexIt = possibility.catalogIndices.begin() + conf*possibility.NumTrueStars();
                        // for (auto iIterator = ijMap.equal_range(catalogIndexIt[*firstNeighborIndexIndex]); iIterator.first != iIterator.second; iIterator.first++) {
                        //     int16_t iCandidate = iIterator.first->second;
                        PairDistanceInvolvingIterator iIterator(ijQuery, ijEnd, catalogIndexIt[*firstNeighborIndexIndex]);
                        while (iIterator.HasValue()) {
                            int16_t iCandidate = *iIterator;
                            ++iIterator;

                            // verify its distance against all other stars in the configuration
                            // TODO: spectral check too

                            // TODO IMPORTANT: (not related to this): We don't really need 4 stars
                            // in range if we only use kvector on the first one. So really, we just
                            // need two stars in range so that the other case works. To summarize,
                            // in the two-star possibility case, we already know what the two
                            // closest are, just use them. In the >2 star possibility case, we only
                            // need the nearest one.
                            for (int whichStar = 0; whichStar < possibility.NumTrueStars(); whichStar++) {
                                float centroidDist = AngleUnit(starSpatials[curCentroid],
                                                               starSpatials[possibility.centroidIndices[whichStar]]);
                                float catalogDist = AngleUnit(catalog[iCandidate].spatial,
                                                              catalog[catalogIndexIt[whichStar]].spatial);

                                // if any star in the configuration is incompatible, the
                                // configuration is incompatible.
                                if (abs(centroidDist - catalogDist) > tolerance) {
                                    goto iContinue;
                                }
                            }

                            // all star distances verified, add the new configuration!
                            // but first, remove it from the old prior: TODO: consider whether this is theoretially justified
                            // Let's be a bit more precise about what a possibility means.
                            // A possibility with a set S of true stars means that the set S is all true, and the set of stars considered so far minus S is all false.
                            // Since the current star "cannot" be false if it's in the same location as a true star, it is in fact correct to remove old configurations in situations like this.
                            oldConfigurationsToRemove.push_back(conf);
                            // add all the stars from the old possibility
                            newPossibility.catalogIndices.insert(newPossibility.catalogIndices.end(),
                                                                 catalogIndexIt, catalogIndexIt + possibility.NumTrueStars());
                            // plus the new star
                            newPossibility.catalogIndices.push_back(iCandidate);

                        iContinue:;
                        }
                    }

                    // remove the old configurations
                    for (int i = oldConfigurationsToRemove.size()-1; i >= 0; i--) {
                        possibility.catalogIndices.erase(possibility.catalogIndices.begin() + oldConfigurationsToRemove[i]*possibility.NumTrueStars(),
                                                         possibility.catalogIndices.begin() + (oldConfigurationsToRemove[i]+1)*possibility.NumTrueStars());
                    }
                }

                if (!newPossibility.catalogIndices.empty()) {
                    prior.push_back(newPossibility);
                }

                break;
            }
            }

            // Possibility of new star being false.
            prior[priorI].probability *= numFalseStars/kUnitSphereArea*probScale;
            // TODO: if there really is a false star, then all possibilities after that point
            // are going to have really small probabilities. Maybe we do need to play around
            // with probscale...
        }

        const int originalNumPossibilities = prior.size();
        const float posteriorSum = BayesPriorTotalProbability(prior, catalog);
        std::cerr << "Unnormalized posterior sum: " << posteriorSum << std::endl;
        probScale *= 1e5/posteriorSum; // TODO: it's probably stable enough, but it's not great:
        std::sort(prior.begin(), prior.end(),
                  // brrt brrt style violation
                  [&catalog](const BayesPossibility &p1, const BayesPossibility &p2) -> bool {
                      // sort in reverse, so highest probability first
                      return p1.TotalProbability(catalog) > p2.TotalProbability(catalog);
                  });
        assert(!prior.empty()); // TODO: can this happen?

        std::cerr << "Debug before trim:" << std::endl;
        //DebugPrintBayesPriorSummary(DebugCalculateBayesPriorSummary(prior, catalog));
        // remove the least likely possibilities, ensuring that the probability of all the
        // removed possibilities doesn't sum to more than admissibleIgnoredProbability
        float discardedSum = 0;
        // this will be zero if there are no configurations, which is good:
        float backProbability = prior.back().TotalProbability(catalog);
        while ((discardedSum + backProbability) / posteriorSum <= admissibleIgnoredProbability) {
            discardedSum += backProbability;
            prior.pop_back();
            assert(!prior.empty());
            backProbability = prior.back().TotalProbability(catalog);
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cerr << "Took " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us" << std::endl;
        std::cerr << "Trimmed from " << originalNumPossibilities << " to " << prior.size() << " possibilities." << std::endl;
#if LOST_DEBUG > 3
        DebugPrintBayesPriorSummary(DebugCalculateBayesPriorSummary(prior, catalog));
#endif

        exploredCentroids.push_back(curCentroid);
    }

    // TODO: it's possible for the mode probability to be below the hard confidence threshold even
    // though we are extremely confident about most stars!
    StarIdentifiers result = Mode(&prior, catalog, hardConfidenceThreshold, stars.size());
    std::cerr << "Idenified " << result.size() << " many stars." << std::endl;
    return result;

    // TODO: use the soft threshold
}

// Bayes Debugging Stuff

// basic summary information about a group of possibilities
class DebugBayesPossibilitiesSummary {
public:
    int numConfigurations = 0;
    // unnormalized
    float totalProbability = 0.0;
};

// the i-th element contains summary for all possibilities with i true stars.
typedef std::vector<DebugBayesPossibilitiesSummary> DebugBayesPriorSummary;

static DebugBayesPriorSummary DebugCalculateBayesPriorSummary(const BayesPrior &prior, const Catalog &catalog) {
    DebugBayesPriorSummary result;
    for (const BayesPossibility &possibility : prior) {
        while ((int)result.size() <= possibility.NumTrueStars()) {
            result.push_back(DebugBayesPossibilitiesSummary());
        }
        result[possibility.NumTrueStars()].numConfigurations += possibility.NumConfigurations(catalog);
        result[possibility.NumTrueStars()].totalProbability += possibility.TotalProbability(catalog);
    }

    return result;
}

// print debug info to stderr
static void DebugPrintBayesPriorSummary(const DebugBayesPriorSummary &summary) {
    float totalProbability = 0.0;
    for (const DebugBayesPossibilitiesSummary &sum : summary) {
        totalProbability += sum.totalProbability;
    }
    for (int i = 0; i < (int)summary.size(); i++) {
        const DebugBayesPossibilitiesSummary &sum = summary[i];
        fprintf(stderr, "%4d true stars: %6d possible configurations, %.5f\n",
                i, sum.numConfigurations, sum.totalProbability/totalProbability);
    }
}

}
