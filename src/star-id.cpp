#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
#include <set>

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
                long numReturnedPairs;
                const int16_t *lowerBoundSearch = vectorDatabase.FindPairsLiberal(
                    lowerBoundRange, upperBoundRange, &numReturnedPairs);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastructure
                for (const int16_t *k = lowerBoundSearch; k < lowerBoundSearch + numReturnedPairs * 2; k++) {
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
        : pairs(NULL), pastTheEnd(NULL) { };

    PairDistanceInvolvingIterator(const int16_t *pairs, long numPairs, int16_t involving)
        : pairs(pairs), pastTheEnd(pairs + numPairs*2), involving(involving) {

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
        return pairs != pastTheEnd;
    }

    // bool operator==(const PairDistanceInvolvingIterator &other) const {
    //     return ()other.pairs == pairs;
    // }

    // bool operator!=(const PairDistanceInvolvingIterator &other) const {
    //     return !(*this == other);
    // }
private:
    const int16_t *pairs;
    const int16_t *pastTheEnd;
    int16_t involving;
    int16_t curValue;

    // like postfix++, except it's a no-op if already on a valid spot.
    void forwardUntilInvolving() {
        while (pairs != pastTheEnd) {
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

// the probability of this match occurring if the stars were uniformly randomly distributed
// throughout the celestial sphere
// float PyramidAnalyticalMatchProbability(float tolerance, int n, )

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

    for (int p = 0; p < stars.size(); p++) {
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
        long ipNum;
        const int16_t *ipPairs = db.FindPairsLiberal(ipDist - tolerance, ipDist + tolerance, &ipNum);
        PairDistanceInvolvingIterator pIterator(ipPairs, ipNum, pyramidIdentifiers[0].catalogIndex);

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
        }
    }
}

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    std::cout << catalog.size() << std::endl;

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL || stars.size() < 4) {
        std::cerr << "Not enough stars, or database missing." << std::endl;
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    // this iteration technique is described in the Pyramid paper. Briefly: i will always be the
    // lowest index, then dj and dk are how many indexes ahead the j-th star is from the i-th, and k-th
    // from the j-th
    long totalIterations = 0;
    for (int dj = 1; dj < (int)stars.size()-1; dj++) {
        for (int dk = 1; dk < (int)stars.size()-dj-1; dk++) {
            for (int dr = 1; dr < (int)stars.size()-dk-dj-1; dr++) {
                for (int i = 0; i < (int)stars.size()-dj-dk; i++) {

                    // identification failure due to cutoff
                    if (++totalIterations > cutoff) {
                        std::cerr << "Cutoff reached." << std::endl;
                        return identified;
                    }

                    int j = i+dj;
                    int k = j+dk;
                    int r = k+dr;

                    // TODO: move this out of the loop?
                    Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
                    Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
                    Vec3 kSpatial = camera.CameraToSpatial(stars[k].position).Normalize();
                    Vec3 rSpatial = camera.CameraToSpatial(stars[r].position).Normalize();

                    float ijDist = AngleUnit(iSpatial, jSpatial);
                    float ikDist = AngleUnit(iSpatial, kSpatial);
                    float irDist = AngleUnit(iSpatial, rSpatial);
                    float jkDist = AngleUnit(jSpatial, kSpatial);
                    float jrDist = AngleUnit(jSpatial, rSpatial);
                    float krDist = AngleUnit(kSpatial, rSpatial);
                    long ijNum, ikNum, irNum; //, jkNum, jrNum, krNum;
                    const int16_t *ijQuery = vectorDatabase.FindPairsLiberal(ijDist - tolerance, ijDist + tolerance, &ijNum);
                    const int16_t *ikQuery = vectorDatabase.FindPairsLiberal(ikDist - tolerance, ikDist + tolerance, &ikNum);
                    const int16_t *irQuery = vectorDatabase.FindPairsLiberal(irDist - tolerance, irDist + tolerance, &irNum);
                    // const int16_t *jkQuery = vectorDatabase.FindPairsLiberal(jkDist - tolerance, jkDist + tolerance, &jkNum);
                    // const int16_t *jrQuery = vectorDatabase.FindPairsLiberal(jrDist - tolerance, jrDist + tolerance, &jrNum);
                    // const int16_t *krQuery = vectorDatabase.FindPairsLiberal(krDist - tolerance, krDist + tolerance, &krNum);

                    PairDistanceInvolvingIterator involvingEnd;
                    std::vector<bool> iCandidates(catalog.size(), false);
                    std::vector<bool> jCandidates(catalog.size(), false);
                    std::vector<bool> kCandidates(catalog.size(), false);
                    std::vector<bool> rCandidates(catalog.size(), false);
                    // TODO: this can be faster, we can combine the initial loop to determine
                    // possible i-s with the loop through them.
                    std::set<int16_t> iAll;
                    for (int p = 0; p < ijNum*2; p++) {
                        iAll.insert(ijQuery[p]);
                    }
                    for (int16_t iCandidate : iAll) {
                        PairDistanceInvolvingIterator jIterator(ijQuery, ijNum, iCandidate);
                        PairDistanceInvolvingIterator kIterator(ikQuery, ikNum, iCandidate);
                        PairDistanceInvolvingIterator rIterator(irQuery, irNum, iCandidate);
                        // TODO: break fast if any of the iterators are empty, if it's any
                        // significant performance improvement.
                        while (jIterator.hasValue()) {
                            while (kIterator.hasValue()) {

                                // small optimization: We can calculate jk before iterating through r, so we will!
                                float jkCandidateDist = AngleUnit(catalog[*jIterator].spatial, catalog[*kIterator].spatial);
                                if (jkCandidateDist < jkDist - tolerance || jkCandidateDist > jkDist + tolerance) {
                                    goto kContinue;
                                }

                                // TODO: if there are no jr matches, there's no reason to
                                // continue iterating through all the other k-s. Possibly
                                // enumarete all r matches, according to ir, before this loop
                                while (rIterator.hasValue()) {
                                    float jrCandidateDist = AngleUnit(catalog[*jIterator].spatial, catalog[*rIterator].spatial);
                                    float krCandidateDist = AngleUnit(catalog[*kIterator].spatial, catalog[*rIterator].spatial);
                                    if (jrCandidateDist < jrDist - tolerance || jrCandidateDist > jrDist + tolerance) {
                                        goto rContinue;
                                    }
                                    if (krCandidateDist < krDist - tolerance || krCandidateDist > krDist + tolerance) {
                                        goto rContinue;
                                    }

                                    // if we made it this far, all 6 angles are confirmed!

                                    // TODO: analytically check confidence

                                    // TODO: rather than immediately proceeding to identify the
                                    // remaining stars, save the identified stars to variables and
                                    // continue the loop to ensure they're unique?
                                    std::cerr << "Surprise muthafuckas" << std::endl;
                                    identified.push_back(StarIdentifier(i, iCandidate));
                                    identified.push_back(StarIdentifier(j, *jIterator));
                                    identified.push_back(StarIdentifier(k, *kIterator));
                                    identified.push_back(StarIdentifier(r, *rIterator));

                                    // TODO: attempt to identify other stars based on the successfully identified pyramid.
                                    PyramidIdentifyRemainingStars(&identified, stars, catalog, vectorDatabase, camera, tolerance);
                                    return identified;

                                rContinue:
                                    ++rIterator;
                                }

                            kContinue:
                                ++kIterator;
                            }

                            ++jIterator;
                        }
                    }                        

                    // int idxI, idxJ, idxK, idxR;;
                    // float confidence = PyramidMatchConfidence(catalog, database,
                    //                                           stars[i], stars[j], stars[k], stars[r],
                    //                                           &idxI, &idxJ, &idxK, &idxR);
                }
            }
        }
    }

    std::cerr << "Tried all pyramids; none matched." << std::endl;
    return identified;
}

StarIdentifiers NonDimStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(TripleDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        return identified;
    }
    TripleDistanceKVectorDatabase vectorDatabase(multiDatabase.SubDatabasePointer(TripleDistanceKVectorDatabase::kMagicValue));

    // lookup to see which image stars have been assigned to which catalog star
    // std::vector<int16_t> identified_fast((int)stars.size(), -1);
    // lookup to see how many times an image star has been identified, -1 if ever misidentified
    // std::vector<int16_t> identified_count((int)stars.size(), 0);
    std::cout << "catalog size: " << catalog.size() << std::endl;
    std::vector<std::map<int16_t, int16_t>> vote((int) catalog.size());
    int16_t num_identified = 0;
    // every possible triangle in the image
    for (int i = 0; i < (int)stars.size(); i++) {  
        for (int j = i + 1; j < (int)stars.size(); j++) {
            for (int k = j + 1; k < (int)stars.size(); k++) {
                /*
                // skip triangle if any of the three stars are misidentified previously
                if (identified_count[i] == -1 || identified_count[j] == -1 || identified_count[k] == -1) {
                    continue;
                }
                */
                int mindex = i;
                int middex = j;
                int maxdex = k;
                // compute target angles
                float smallAngle;
                float midAngle;
                float largeAngle;
                FocalPlaneAngles(stars, camera, smallAngle, midAngle, largeAngle, mindex, middex, maxdex, i, j, k);
                if (std::isnan(smallAngle) || std::isnan(midAngle) || std::isnan(largeAngle)) {
                    continue;
                }
                // range of query
                float lowerBoundRange = smallAngle - tolerance;
                float upperBoundRange = smallAngle + tolerance;
                long numReturnedTriples;
                const int16_t *lowerBoundSearch = vectorDatabase.FindTriplesLiberal(
                    lowerBoundRange, upperBoundRange, &numReturnedTriples);
                if (numReturnedTriples < 1) {
                    continue;
                }
                // matched triple
                const int16_t* mt = NULL;
                bool unique = true;
                /*
                 * basically if we make a query on small angle, all returned catalog triangles have small angle
                 * within tolerance and we want to compare each resulting triangle to affirm the large angle is 
                 * also within tolerance but if there is more than one than we skip this query
                */
                for (const int16_t *l = lowerBoundSearch; l < lowerBoundSearch + numReturnedTriples * 3; l += 3) {
                    int actualMindex;
                    int actualMiddex;
                    int actualMaxdex;
                    float actualSmallAngle;
                    float actualMidAngle;
                    float actualLargeAngle;
                    InnerAngles(catalog, actualSmallAngle, actualMidAngle, actualLargeAngle, actualMindex, actualMiddex, actualMaxdex, *l, *(l+1), *(l+2));
                    if (actualLargeAngle > largeAngle - tolerance && actualLargeAngle < largeAngle + tolerance && actualMidAngle > midAngle - tolerance &&
                    actualMidAngle < midAngle + tolerance) {
                        if (mt != NULL) {
                            unique = false;
                        }
                        mt = l;
                    }
                }
                // if there was no unique match then find new target triangle
                if (mt == NULL || !unique) {
                    continue;
                }
                
                // figure out which index is which (compare angles small medium large)
                int actualMaxdex;
                int actualMiddex;
                int actualMindex;
                // garbage chute (AKA dont care)
                float gc;
                InnerAngles(catalog, gc, gc, gc, actualMindex, actualMiddex, actualMaxdex, *mt, *(mt+1), *(mt+2));
                if (!vote[actualMindex].count(mindex)) {
                    vote[actualMindex][mindex] = 1;
                } else {
                    vote[actualMindex][mindex]++;
                }
                if (!vote[actualMiddex].count(middex)) {
                    vote[actualMiddex][middex] = 1;
                } else {
                    vote[actualMiddex][middex]++;
                }
                if (!vote[actualMaxdex].count(maxdex)) {
                    vote[actualMaxdex][maxdex] = 1;
                } else {
                    vote[actualMaxdex][maxdex]++;
                }
                /*
                // confirm stars do not get misidentified and identify these three stars 
                if (identified_fast[mindex] != -1 && identified_fast[mindex] != actualMindex) {
                    identified_count[mindex] = -1;
                }
                if (identified_fast[middex] != -1 && identified_fast[middex] != actualMiddex) {
                    identified_count[middex] = -1;
                }
                if (identified_fast[maxdex] != -1 && identified_fast[maxdex] != actualMaxdex) {
                    identified_count[maxdex] = -1;
                }
                // if each of the three stars has not been misidentified
                if (identified_count[mindex] != -1 && identified_count[middex] != -1 && identified_count[maxdex] != -1) {
                    // increment identification count, set their identification (may already be set)
                    identified_count[mindex]++;
                    identified_count[middex]++;
                    identified_count[maxdex]++;
                    if (identified_count[mindex] == 1 + num_verify) {
                        num_identified++;
                    }
                    if (identified_count[middex] == 1 + num_verify) {
                        num_identified++;
                    }
                    if (identified_count[maxdex] == 1 + num_verify) {
                        num_identified++;
                    }
                    if (num_identified >= 20) {
                        for (int i = 0; i < (int)stars.size(); i++) {
                            if (identified_count[i] >= 1 + num_verify) {
                                StarIdentifier newStar(i, identified_fast[i]);
                                identified.push_back(newStar);
                            }
                        }
                        return identified;
                    }
                    identified_fast[mindex] = actualMindex;
                    identified_fast[middex] = actualMiddex;
                    identified_fast[maxdex] = actualMaxdex;
                }
                */
            }
        }
    }
    std::cout << "Made it here" << std::endl;
    for (int i = 0; i < (int) catalog.size(); i++) {
        auto pr = std::max_element(std::begin(vote[i]), std::end(vote[i]),[] (const std::pair<int16_t,int16_t>& a, const std::pair<int16_t,int16_t>& b)->bool{ return a.second < b.second; } );
        int16_t maxval = pr->second;
        int16_t maxkey = pr->first;
        if (maxval > 0 + num_verify) {
            StarIdentifier newStar(maxkey, i);
            identified.push_back(newStar);
        }
    }
    // finalize identification of each image star that has been identified at least once plus number of verifications
    // without multiple identifications to different stars (misidentification)
    /*
    for (int i = 0; i < (int)stars.size(); i++) {
        if (identified_count[i] >= 1 + num_verify) {
            StarIdentifier newStar(i, identified_fast[i]);
            identified.push_back(newStar);
        }
    }
    */
    return identified;
}

}
