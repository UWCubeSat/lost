#include <limits>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "star-id.hpp"
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

/// unidentified centroid used in IdentifyRemainingStarsPairDistance
/// The "angles" through here are "triangular angles". A triangular angle is a 2D angle in a triangle formed by three centroids.
class IRUnidentifiedCentroid {
public:
    IRUnidentifiedCentroid(const Star &star)
        : bestAngleFrom90(std::numeric_limits<float>::max()), // should be infinity
          bestStar1(0,0), bestStar2(0,0),
          star(star) { }

    float bestAngleFrom90; /// For the pair of other centroids forming the triangular angle closest to 90 degrees, how far from 90 degrees it is (in radians)
    StarIdentifier bestStar1; /// One star corresponding to bestAngleFrom90
    StarIdentifier bestStar2; /// The other star corresponding to bestAngleFrom90
    int16_t index; /// Index into list of all centroids
    const Star &star;

private:
    // possible improvement: Use a tree map here to allow binary search
    std::vector<std::pair<float, StarIdentifier>> identifiedStarsInRange;

private:
    float VerticalAnglesToAngleFrom90(float v1, float v2) {
        return abs(FloatModulo(v1-v2, M_PI) - M_PI_2);
    }

public:
    /**
     * When a centroid within range of this centroid is identified, call this function. This
     * function does /not/ check whether the centroid is within range.
     */
    void AddIdentifiedStar(const StarIdentifier &starId, const Stars &stars) {
        const Star &otherStar = stars[starId.starIndex];
        Vec2 positionDifference = otherStar.position - star.position;
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
};

std::vector<IRUnidentifiedCentroid *> FindUnidentifiedCentroidsInRange(std::vector<IRUnidentifiedCentroid> *centroids, const Star &star, float minDistance, float maxDistance) {
    std::vector<IRUnidentifiedCentroid *> result;

    // Find the first centroid that is within range of the given centroid.
    auto firstInRange = std::lower_bound(centroids->begin(), centroids->end(), star.position.x - maxDistance,
        [](const IRUnidentifiedCentroid &centroid, float x) {
            return centroid.star.position.x < x;
        });

    // Find the first centroid that is not within range of the given centroid.
    auto firstNotInRange = std::lower_bound(firstInRange, centroids->end(), star.position.x + maxDistance,
        [](const IRUnidentifiedCentroid &centroid, float x) {
            return centroid.star.position.x <= x;
        });

    // Copy the pointers to the stars into the result vector.
    for (auto it = firstInRange; it != firstNotInRange; ++it) {
        float distance = Distance(star.position, it->star.position);
        if (distance >= minDistance && distance <= maxDistance) {
            result.push_back(&*it);
        }
    }

    return result;
}

/**
 * Given a newly identified centroid and a list of unidentifed centroids, update all the
 * unidentified centroids which are within range.
 */
void AddToAllUnidentifiedCentroids(const StarIdentifier &starId, const Stars &stars,
                                   std::vector<IRUnidentifiedCentroid> *unidentifiedCentroids,
                                   float minDistance, float maxDistance,
                                   const Camera &camera) {
    for (IRUnidentifiedCentroid *centroid : FindUnidentifiedCentroidsInRange(unidentifiedCentroids, stars[starId.starIndex], minDistance, maxDistance)) {
        centroid->AddIdentifiedStar(starId, stars);
    }
}

/**
 * Given some identified stars, attempt to identify the rest.
 *
 * Requires a pair distance database to be present. Iterates through the unidentified centroids in
 * an intelligent order, identifying them one by one.
 */
int IdentifyRemainingStarsPairDistance(StarIdentifiers *identifiers,
                                        const Stars &stars,
                                        const Catalog &catalog,
                                        const PairDistanceKVectorDatabase &db,
                                        const Camera &camera,
                                        float tolerance) {
    // initialize all unidentified centroids
    std::vector<IRUnidentifiedCentroid> unidentifiedCentroids;
    for (const Star &star : stars) {
        unidentifiedCentroids.push_back(IRUnidentifiedCentroid(star));
    }

    // sort unidentified centroids by x coordinate
    std::sort(unidentifiedCentroids.begin(), unidentifiedCentroids.end(),
        [](const IRUnidentifiedCentroid &a, const IRUnidentifiedCentroid &b) {
            return a.star.position.x < b.star.position.x;
        });

    // for each identified star, add it to the list of identified stars for each unidentified centroid within range
    for (const auto &starId : *identifiers) {
        for (auto it = unidentifiedCentroids.begin(); it != unidentifiedCentroids.end();) {
            if (starId.starIndex == it->index) {
                unidentifiedCentroids.erase(it);
            } else {
                AddToAllUnidentifiedCentroids(starId, stars, &unidentifiedCentroids, db.MinDistance(), db.MaxDistance(), camera);
            }
        }
    }

    int numExtraIdentifiedStars = 0;

    // keep getting the best unidentified centroid and identifying it
    while (!unidentifiedCentroids.empty()) {
        // find the best unidentified centroid
        auto bestUnidentifiedCentroidIt = std::min_element(unidentifiedCentroids.begin(), unidentifiedCentroids.end(),
            [](const IRUnidentifiedCentroid &a, const IRUnidentifiedCentroid &b) {
                return a.bestAngleFrom90 < b.bestAngleFrom90;
            });
        const IRUnidentifiedCentroid &bestUnidentifiedCentroid = *bestUnidentifiedCentroidIt;

        // Project best stars to 3d, find angle between them and current unidentified centroid
        Vec3 unidentifiedSpatial = camera.CameraToSpatial(bestUnidentifiedCentroid.star.position);
        float d1 = Angle(camera.CameraToSpatial(stars[bestUnidentifiedCentroid.bestStar1.starIndex].position), unidentifiedSpatial);
        float d2 = Angle(camera.CameraToSpatial(stars[bestUnidentifiedCentroid.bestStar2.starIndex].position), unidentifiedSpatial);

        const int16_t *query1End, *query2End;
        const int16_t *query1 = db.FindPairsLiberal(d1-tolerance, d1+tolerance, &query1End);
        const int16_t *query2 = db.FindPairsLiberal(d2-tolerance, d2+tolerance, &query2End);

        // Use PairDistanceInvolvingIterator to find catalog candidates for the unidentified centroid from both sides.
        PairDistanceInvolvingIterator it1(query1, query1End, bestUnidentifiedCentroid.bestStar1.catalogIndex);
        PairDistanceInvolvingIterator it2(query2, query2End, bestUnidentifiedCentroid.bestStar2.catalogIndex);

        // find all the catalog stars that are in both annuli
        std::vector<int16_t> candidates1 = ConsumeInvolvingIterator(it1);
        std::vector<int16_t> candidates2 = ConsumeInvolvingIterator(it2);
        std::sort(candidates1.begin(), candidates1.end());
        std::sort(candidates2.begin(), candidates2.end());
        std::vector<int16_t> candidatesIntersection;
        std::set_intersection(candidates1.begin(), candidates1.end(), candidates2.begin(), candidates2.end(), std::back_inserter(candidatesIntersection));

        if (candidatesIntersection.size() != 1) { // if there is not exactly one candidate, we can't identify the star. Just remove it from the list.
            if (candidatesIntersection.size() > 1) {
                std::cerr << "WARNING: Multiple catalog stars matched during identify remaining stars. This should be rare." << std::endl;
            }
        } else {
            // identify the centroid
            identifiers->emplace_back(bestUnidentifiedCentroid.index, candidatesIntersection[0]);

            // update nearby unidentified centroids with the new identified star
            AddToAllUnidentifiedCentroids(identifiers->back(), stars, &unidentifiedCentroids, db.MinDistance(), db.MaxDistance(), camera);

            ++numExtraIdentifiedStars;
        }
        // remove the centroid from the list of unidentified centroids
        unidentifiedCentroids.erase(bestUnidentifiedCentroidIt);
    }

    return numExtraIdentifiedStars;
}

/// After some stars have been identified, try to idenify the rest using a faster algorithm.
void PyramidIdentifyRemainingStars(StarIdentifiers *identifiers,
                                   const Stars &stars,
                                   const Catalog &catalog,
                                   const PairDistanceKVectorDatabase &db,
                                   const Camera &camera,
                                   float tolerance) {

    assert(identifiers->size() == 4);
    StarIdentifiers pyramidIdentifiers = *identifiers; // copy with only the pyramid's high confidence stars
    Vec3 pyramidActualSpatials[4];
    for (int l = 0; l < 4; l++) {
        pyramidActualSpatials[l] = camera.CameraToSpatial(stars[pyramidIdentifiers[l].starIndex].position).Normalize();
    }

    for (int p = 0; p < (int)stars.size(); p++) {
        // ensure this star isn't in the pyramid
        bool pInPyramid = false;
        for (const StarIdentifier &id : pyramidIdentifiers) {
            if (id.starIndex == p) {
                pInPyramid = true;
                break;
            }
        }
        if (pInPyramid) {
            continue;
        }

        Vec3 pSpatial = camera.CameraToSpatial(stars[p].position).Normalize();
        float ipDist = AngleUnit(pyramidActualSpatials[0], pSpatial);
        const int16_t *ipEnd;
        const int16_t *ipPairs = db.FindPairsLiberal(ipDist - tolerance, ipDist + tolerance, &ipEnd);
        PairDistanceInvolvingIterator pIterator(ipPairs, ipEnd, pyramidIdentifiers[0].catalogIndex);

        std::vector<int16_t> pCandidates; // collect them all in the loop, at the end only identify
                                          // the star if unique
        while (pIterator.HasValue()) {
            bool ok = true;
            for (int l = 1; l < 4; l++) {
                float actualDist = AngleUnit(pSpatial, pyramidActualSpatials[l]);
                float expectedDist = AngleUnit(catalog[*pIterator].spatial,
                                               catalog[pyramidIdentifiers[l].catalogIndex].spatial);
                if (actualDist < expectedDist - tolerance || actualDist > expectedDist + tolerance) {
                    ok = false;
                }
            }
            if (ok) {
                pCandidates.push_back(*pIterator);
            }
            ++pIterator;
        }

        if (pCandidates.size() == 1) {
            identifiers->push_back(StarIdentifier(p, pCandidates[0]));
        }
        if (pCandidates.size() > 1) {
            std::cerr << "duplicate other star??" << std::endl;
            for (int16_t c : pCandidates) {
                std::cerr << catalog[c].name << std::endl;
            }
        }
    }
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


                    int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
                    std::vector<bool> iSeen(catalog.size(), false);
                    for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd; iCandidateQuery++) {
                        int iCandidate = *iCandidateQuery;
                        if (iSeen[iCandidate]) {
                            continue;
                        }
                        iSeen[iCandidate] = true;

                        const Vec3 &iCandidateSpatial = catalog[iCandidate].spatial;

                        // TODO: caching these iterator results into vectors can improve
                        // performance, but at the cost of memory. It would be best to put some kind
                        // of guarantee on the memory usage, and then switch to using the iterator
                        // without caching if that memory limit is exceeded.
                        PairDistanceInvolvingIterator jIterator(ijQuery, ijEnd, iCandidate);
                        PairDistanceInvolvingIterator kIterator(ikQuery, ikEnd, iCandidate);
                        PairDistanceInvolvingIterator rIterator(irQuery, irEnd, iCandidate);
                        std::vector<int16_t> jCandidates = ConsumeInvolvingIterator(jIterator);
                        std::vector<int16_t> kCandidates = ConsumeInvolvingIterator(kIterator);
                        std::vector<int16_t> rCandidates = ConsumeInvolvingIterator(rIterator);
                        // TODO: break fast if any of the iterators are empty, if it's any
                        // significant performance improvement.
                        for (int16_t jCandidate : jCandidates) {
                            const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;
                            Vec3 ijCandidateCross = iCandidateSpatial.CrossProduct(jCandidateSpatial);

                            for (int16_t kCandidate : kCandidates) {
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
                                for (int16_t rCandidate : rCandidates) {
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
                    }

                    if (iMatch != -1) {
                        printf("Matched unique pyramid!\nExpected mismatches: %e\n", expectedMismatches);
                        identified.push_back(StarIdentifier(i, iMatch));
                        identified.push_back(StarIdentifier(j, jMatch));
                        identified.push_back(StarIdentifier(k, kMatch));
                        identified.push_back(StarIdentifier(r, rMatch));

                        PyramidIdentifyRemainingStars(&identified, stars, catalog, vectorDatabase, camera, tolerance);
                        printf("Identified an additional %d stars\n", (int)identified.size() - 4);

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

}
