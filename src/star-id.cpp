#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

#include "star-id.hpp"
#include "databases.hpp"
#include "attitude-utils.hpp"

namespace lost {

StarIdentifiers DummyStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers result;

    for (int i = 0; i < (int)stars.size(); i++) {
        if (rand() > RAND_MAX/5) {
            result.push_back(StarIdentifier(i, rand() % 2000));
        }
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

    /**
     * Strategies:
     * 
     * 1. For each star, enumerate all stars which have the same combination of distances to some
     *  other stars, getting down to a hopefully small (<10) list of candidates for each star, then
     *  do a quad-nested loop to correlate them.
     *
     * 2. Loop through all possible stars in the catalog for star i. Then look at edge ij, using
     * this to select possible j-th stars. If ever there is not a possible j-th star, continue the
     * i-loop. When a possible ij combination is found, loop through k stars according to ik. IF
     * none are found, continue the outer i loop. If some are found, check jk for each one. For each possible ijk triangle, 
     */

class PairDistanceInvolvingIterator {
public:
    // unqualified constructor makes a "past-the-end" iterator
    PairDistanceInvolvingIterator()
        : pairs(NULL), end(NULL) { };

    PairDistanceInvolvingIterator(const int16_t *pairs, const int16_t *end, int16_t involving)
        : pairs(pairs), end(end), involving(involving) {

        assert((end-pairs)%2 == 0);
        forwardUntilInvolving();
    };

    // PairDistanceInvolvingIterator operator++() {
    //     PairDistanceInvolvingIterator result(*this);
    //     ++(*this);
    //     return result;
    // }

    PairDistanceInvolvingIterator &operator++() {
        assert(hasValue());
        pairs += 2;
        forwardUntilInvolving();
        return *this;
    }

    int16_t operator*() const {
        return curValue;
    }

    bool hasValue() {
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

    // like postfix++, except it's a no-op if already on a valid spot.
    void forwardUntilInvolving() {
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
        while (pIterator.hasValue()) {
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

                    assert(i!=j && j!=k && k!=r && i!=k && i!=r && j!=r);

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
                    bool spectralTorch = iSpatial.crossProduct(jSpatial)*kSpatial > 0;

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
                        std::vector<int16_t> jCandidates;
                        std::vector<int16_t> kCandidates;
                        std::vector<int16_t> rCandidates;
                        while (jIterator.hasValue()) {
                            jCandidates.push_back(*jIterator);
                            ++jIterator;
                        }
                        while (kIterator.hasValue()) {
                            kCandidates.push_back(*kIterator);
                            ++kIterator;
                        }
                        while (rIterator.hasValue()) {
                            rCandidates.push_back(*rIterator);
                            ++rIterator;
                        }
                        // TODO: break fast if any of the iterators are empty, if it's any
                        // significant performance improvement.
                        for (int16_t jCandidate : jCandidates) {
                            const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;
                            Vec3 ijCandidateCross = iCandidateSpatial.crossProduct(jCandidateSpatial);

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
                     const int16_t *queryBegin, const int16_t *queryEnd)
        : probability(probability) {
        centroidIndices.push_back(centroidIndex1);
        centroidIndices.push_back(centroidIndex2);
        catalogIndicesBegin = queryBegin;
        catalogIndicesEnd = queryEnd;
    }
    BayesPossibility() = default;

    // the number of subconfigurations stored in this possibility (subconfigurations are really just
    // more possibilities)
    int NumConfigurations(const Catalog &catalog) const {
        if (catalogIndicesBegin != NULL) {
            assert(centroidIndices.size() == 2);
            // 2* to account for flipping
            return 2*(catalogIndicesEnd - catalogIndicesBegin) / centroidIndices.size();
        } else if (catalogIndices.size() == 1) {
            return catalog.size();
        } else if (catalogIndices.size() == 0) {
            return 1;
        } else {
            return catalogIndices.size() / centroidIndices.size();
        }
    }

    int NumTrueStars() const {
        return centroidIndices.size();
    }

    float TotalProbability(const Catalog &catalog) const {
        return probability * NumConfigurations(catalog);
    }

    void MultiplyProbabilityBy(float rescale) {
        probability *= rescale;
    }

    // convert this possibility to an identifier if appropriate
    StarIdentifiers ToStarIdentifiers() const {
        assert(centroidIndices.size() > 2);
        assert(centroidIndices.size() == catalogIndices.size());
        StarIdentifiers result;
        for (int i = 0; i < centroidIndices.size(); i++) {
            result.push_back(StarIdentifier(centroidIndices[i], catalogIndices[i]));
        }
        return result;
    };

    float probability; // Probability of each possible configuration inside
    std::vector<int16_t> centroidIndices;
    const int16_t *catalogIndicesBegin = NULL;
    const int16_t *catalogIndicesEnd;
    std::vector<int16_t> catalogIndices;
};

// can iterate over either manually passed in int16_t*s, or vector iterators. When iterating over
// int16_t*s, also goes over them in reverse
class BayesPossibilityIterator {
public:
    BayesPossibilityIterator(BayesPossibility possibility)
        : queryStart(possibility.catalogIndicesBegin), queryCur(possibility.catalogIndicesBegin), queryEnd(possibility.catalogIndicesEnd),
          vecCur(possibility.catalogIndices.cbegin()), vecEnd(possibility.catalogIndices.cend())
        { };

    bool hasValue() const {
        if (queryCur != NULL) {
            return !(queryCur == queryStart && reverse == true);
        } else {
            return vecCur != vecEnd;
        }
    }
    void operator++() {
        assert(hasValue());
        if (queryCur != NULL) {
            if (!reverse) {
                queryCur++;
                if (queryCur == queryEnd) {
                    queryCur--;
                    reverse = true;
                }
            } else {
                queryCur--;
            }
        } else {
            vecCur++;
        }
    };
    int16_t operator*() {
        assert(hasValue());
        if (queryCur != NULL) {
            return *queryCur;
        } else {
            return *vecCur;
        }
    };

private:
    const int16_t *queryStart;
    const int16_t *queryCur = NULL;
    const int16_t *queryEnd;
    bool reverse = false;
    std::vector<int16_t>::const_iterator vecCur;
    std::vector<int16_t>::const_iterator vecEnd;
};

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

// return the mode of the distribution. If nothing promising, return empty and negative probability
StarIdentifiers Mode(const BayesPrior &prior, const Catalog &catalog, float *modeProbability) {
    float sum = 0.0;
    float max = -1.0;
    BayesPossibility maxPossibility;
    for (const BayesPossibility &possibility : prior) {
        sum += possibility.TotalProbability(catalog);
        // TODO: check that numconfigurations always storable in 16 bits? Hopefully should be, but
        // maybe just assert?
        int numConfigurations = possibility.NumConfigurations(catalog);
        if (numConfigurations == 1 && possibility.probability > max) {
            max = possibility.probability;
            maxPossibility = possibility;
        }
    }
    if (max < 0) {
        *modeProbability = -1.0;
        return StarIdentifiers();
    } else {
        *modeProbability = max/sum;
        return maxPossibility.ToStarIdentifiers();
    }
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

StarIdentifiers BayesianStarIdAlgorithm::Go(unsigned char *database, const Stars &stars,
                                            const Catalog &catalog, const Camera &camera) const {
    const float kUnitSphereArea = 4*M_PI;
    // for each new centroid, we consider information from up to kMaxNearbyCentroids many centroids
    const int kIdealNumNeighbors = 4;
    const int kMinNumNeighbors = 2;
    // if there are fewer than this many possibilities for a star, check them all individually
    // rather than using an involving iterator to narrow them down.
    const int kInvolvingIteratorCutoff = 30;

    // load database
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        std::cerr << "Database missing" << std::endl;
        return StarIdentifiers();
    }
    PairDistanceKVectorDatabase kvector(databaseBuffer);

    std::vector<Vec3> starSpatials;
    for (const Star &star : stars) {
        starSpatials.push_back(camera.CameraToSpatial(star.position).Normalize());
    }

    // we only keep track of un-normalized probabilities (do not sum to one). So, we can multiply
    // all of them by an arbitrary constant without functionally changing anything. This variable is
    // picked just to keep probabilities in a reasonable range. Towards this end, we set constant
    // approximately so that multiplying by 
    const float probScale = tolerance*tolerance*3*M_PI;
    const float falseStarProbability = (float)numFalseStars/kUnitSphereArea;

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

        for (int centroidOffset = 0; centroidOffset < stars.size(); centroidOffset++) {
            const int curCentroid = (firstCentroid + centroidOffset) % stars.size();

            // ensure that there are enough neighbors
            int numNeighbors = 0;
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
                    if (numNeighbors >= minNumNeighbors) {
                        break;
                    }
                }
            }
            if (alreadyExplored || numNeighbors < minNumNeighbors) {
                continue;
            }
            // There are enough neighbors!
            searchMode = CentroidSearchMode::Found;
            exploredCentroids.push_back(curCentroid);

                    // find the closest numNeighbors neighbors. Closest will make the rest of this
                    // faster and tend to be well distributed around different sides of the new
                    // centroid, which keeps the annulus intersection area small.
                    // brrt brrt style guideline violation: no lambdas
                    std::nth_element(exploredCentroids.begin(), exploredCentroids.begin() + numNeighbors - 1, exploredCentroids.end(),
                                     [curCentroid, &stars](int s1, int s2) -> bool {
                                         const float s1Dist = (stars[curCentroid].position - stars[s1].position).MagnitudeSq();
                                         const float s2Dist = (stars[curCentroid].position - stars[s2].position).MagnitudeSq();
                                         return s1Dist < s2Dist;
                                     });
                    // nth_element sorts the array up to the given mid iterator.

            // posterior states where the new star is true. We add false states at the end when
            // appropriate.
            BayesPrior posterior;

            for (const BayesPossibility &possibility : prior) {
                // // a more theoretically accurate way would not include possibilities where the current
                // // star must be true. Ie, any configuration that indicates a true star at the current
                // // location certainly shouldn't be included. TODO wouldn't be hard to improve
                // priorSum += possibility.TotalProbability(catalog);

                switch (possibility.centroidIndices.size()) {
                case 0: {
                    float trueProb = possibility.probability/kUnitSphereArea * probScale;
                    posterior.push_back(BayesPossibility(trueProb, curCentroid));

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
                    if (begin != end) {
                        // area of spherical annulus with radius sigma
                        const float annulusArea = 2*M_PI*tolerance*sin(distance);
                        const float trueProb = possibility.probability/annulusArea * probScale;
                        posterior.push_back(
                            BayesPossibility(
                                trueProb, possibility.centroidIndices[0], curCentroid, begin, end));
                    }
                    break;
                }

                default: { // >= 2
                    assert(numNeighbors >= 2);
                    // TODO: handle situation  where the first one is too close (<MinDistance())

                    // TODO: what to do if this property isn't satisfied?

                    /*
                     * When exactly two neighbors, 
                     */
                    BayesPossibilityIterator catalogIndexIt(possibility);
                    int16_t firstCentroid = -1;
                    while (catalogIndexIt.hasValue()) {
                        for (int16_t centroidIndex : possibility.centroidIndices) {
                            float angle = AngleUnit(starSpatials[centroidIndex], starSpatials[curCentroid]);
                            if (angle > kvector.MinDistance()+tolerance && angle < kvector.MaxDistance()-tolerance) {
                                if (firstCentroid < 0) {
                                    firstCentroid = centroidIndex;
                                } else {
                                    // we have two centroids. Do the thing.
                                    // TODO: maybe be pickier about which centroid to pick when there are more than two options
                                    const int16_t *candidatesEnd;
                                    const int16_t *candidates = kvector.FindPairsLiberal(angle-tolerance, angle+tolerance, &candidatesEnd);
                                    PairDistanceInvolvingIterator subCandidates(candidates, candidatesEnd, firstCentroid);
                                }
                            }
                            ++catalogIndexIt;
                        }
                    }
                }
                }

                // Possibility of new star being false.
                BayesPossibility falsePossibility = possibility;
                falsePossibility.MultiplyProbabilityBy(numFalseStars/kUnitSphereArea*probScale);
                posterior.push_back(falsePossibility);
                // TODO: if there really is a false star, then all possibilities after that point
                // are going to have really small probabilities. Maybe we do need to play around
                // with probscale...
            }

            const int originalNumPossibilities = posterior.size();
            const float posteriorSum = BayesPriorTotalProbability(posterior, catalog);
            std::cerr << "Unnormalized posterior sum: " << posteriorSum << std::endl;
            std::sort(posterior.begin(), posterior.end(),
                      // brrt brrt style violation
                      [&catalog](const BayesPossibility &p1, const BayesPossibility &p2) -> bool {
                          // sort in reverse, so highest probability first
                          return p1.TotalProbability(catalog) > p2.TotalProbability(catalog);
                      });
            assert(!posterior.empty()); // TODO: can this happen?
            float discardedSum = 0;
            float backProbability = posterior.back().TotalProbability(catalog);
            while ((discardedSum + backProbability) / posteriorSum <= admissibleIgnoredProbability) {
                discardedSum += posterior.back().TotalProbability(catalog);
                posterior.pop_back();
                assert(!posterior.empty());
                backProbability = posterior.back().TotalProbability(catalog);
            }
            std::cerr << "Trimmed from " << originalNumPossibilities << " to " << posterior.size() << " possibilities." << std::endl;
            // // decide whether we need to consider the current star being false
            // if (posteriorFalseSum/posteriorSum > admissibleIgnoredProbability) {
            //     std::cerr << "False star considered for " << exploredCentroids.size() << std::endl;
            //     // and if so, add all possibilities from the prior back in, properly adjusting their
            //     // probabilities to take into account that the new star is false.
            //     for (const BayesPossibility &priorPossibility : prior) {
            //         BayesPossibility posteriorPossibility = priorPossibility;
            //         posteriorPossibility.MultiplyProbabilityBy(numFalseStars / kUnitSphereArea * probScale);
            //         posterior.push_back(std::move(posteriorPossibility));
            //     }
            // }

            // the ol' switcheroo
            prior = std::move(posterior);
        }
    }

    // TODO: it's possible for the mode probability to be below the hard confidence threshold even
    // though we are extremely confident about most stars!
    float modeProbability;
    StarIdentifiers result = Mode(prior, catalog, &modeProbability);
    if (modeProbability >= hardConfidenceThreshold) {
        return result;
    } else {
        std::cerr << "Didn't meet confidence threshold: " << modeProbability << std::endl;
        return StarIdentifiers();
    }

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
        while (result.size() <= possibility.NumTrueStars()) {
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
