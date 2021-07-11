#include <stdlib.h>
#include <math.h>
#include <assert.h>

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
    PairDistanceKVectorDatabase vectorDatabase(multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue));

    for (int i = 0; i < (int)stars.size(); i++) {  
        std::vector<int16_t> votes(catalog.size(), 0);
        Vec3 iSpatial = camera.CameraToSpatial({ stars[i].x, stars[i].y }).Normalize();
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                // TODO: find a faster way to do this:
                std::vector<bool> votedInPair(catalog.size(), false);
                Vec3 jSpatial = camera.CameraToSpatial({ stars[j].x, stars[j].y }).Normalize();
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
            Vec3 firstSpatial = camera.CameraToSpatial({firstIdentified.x, firstIdentified.y});
            Vec3 secondSpatial = camera.CameraToSpatial({secondIdentified.x, secondIdentified.y});
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

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {
    // TODO
        StarIdentifiers identified;
        return identified;
    ;
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
    int16_t identified_fast[(int)stars.size()] {-1};
    // lookup to see how many times an image star has been identified, -1 if ever misidentified
    int16_t identified_count[(int)stars.size()] {0};

    // every possible triangle in the image
    for (int i = 0; i < (int)stars.size(); i++) {  
        for (int j = i + 1; j < (int)stars.size(); j++) {
            for (int k = j + 1; k < (int)stars.size(); k++) {
                // skip triangle if any of the three stars are misidentified previously
                if (identified_count[i] == -1 || identified_count[j] == -1 || identified_count[k] == -1) {
                    continue;
                }
                // compute target small angle
                int mindex = i;
                int middex = j;
                int maxdex = k;
                float smallAngle = minFocalPlaneAngle(stars, mindex, i, j, k);
                // compute target large angle
                float largeAngle = maxFocalPlaneAngle(stars, maxdex, i, j, k);
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
                 * within toleranceand we want to compare each resulting triangle to affirm the large angle is 
                 * also within tolerance but if there is more than one than we skip this query
                */
                for (const int16_t *l = lowerBoundSearch; l < lowerBoundSearch + numReturnedTriples * 3; l += 3) {
                    int actualMaxdex;
                    float actualLargeAngle = maxInnerAngle(catalog, actualMaxdex, *l, *(l+1), *(l+2));
                    if (actualLargeAngle > largeAngle - tolerance && actualLargeAngle < largeAngle + tolerance) {
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
                
                // figure out which one is which (compare angles small medium large)
                if (i != mindex && i != maxdex) {
                    middex = i;
                } else if (j != mindex && j != maxdex) {
                    middex = j;
                } else {
                    middex = k;
                }
                int actualMaxdex;
                int actualMiddex;
                int actualMindex;
                maxInnerAngle(catalog, actualMaxdex, *mt, *(mt+1), *(mt+2));
                minInnerAngle(catalog, actualMindex, *mt, *(mt+1), *(mt+2));
                if (*mt != actualMindex && *mt != actualMaxdex) {
                    actualMiddex = *mt;
                } else if (*(mt+1) != actualMindex && *(mt+1) != actualMaxdex) {
                    actualMiddex = *(mt+1);
                } else {
                    actualMiddex = *(mt+2);
                }
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
                    identified_fast[mindex] = actualMindex;
                    identified_fast[middex] = actualMiddex;
                    identified_fast[maxdex] = actualMaxdex;
                }
            }
        }
    }
    // finalize identification of each image star that has been identified at least thrice
    // without multiple identifications to different stars (misidentification)
    for (int i = 0; i < (int)stars.size(); i++) {
        if (identified_count[i] >= 3) {
            StarIdentifier newStar(i, identified_fast[i]);
            identified.push_back(newStar);
        }
    }
    return identified;
}

}
