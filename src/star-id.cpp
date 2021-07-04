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

    // constant lookup to see if one of the stars in a triangle has already been identified
    char identified_fast[(int)stars.size()] {0};

    // every possible triangle in the image
    for (int i = 0; i < (int)stars.size(); i++) {  
        for (int j = i + 1; j < (int)stars.size(); j++) {
            for (int k = j + 1; k < (int)stars.size(); k++) {
                if (identified_fast[i] || identified_fast[j] || identified_fast[k]) {
                    continue;
                }
                // compute target small angle
                float smallAngle = std::min(std::min(Angle(catalog[j].spatial-catalog[i].spatial, 
            catalog[k].spatial-catalog[i].spatial), Angle(catalog[j].spatial-catalog[k].spatial, 
            catalog[j].spatial-catalog[i].spatial)), Angle(catalog[k].spatial-catalog[i].spatial, 
            catalog[k].spatial-catalog[j].spatial));
                // compute target large angle
                float largeAngle = std::max(std::max(Angle(catalog[j].spatial-catalog[i].spatial, 
            catalog[k].spatial-catalog[i].spatial), Angle(catalog[j].spatial-catalog[k].spatial, 
            catalog[j].spatial-catalog[i].spatial)), Angle(catalog[k].spatial-catalog[i].spatial, 
            catalog[k].spatial-catalog[j].spatial));
                // range of query
                float lowerBoundRange = smallAngle - tolerance;
                float upperBoundRange = smallAngle + tolerance;
                long numReturnedTriples;
                const int16_t *lowerBoundSearch = vectorDatabase.FindTriplesLiberal(
                    lowerBoundRange, upperBoundRange, &numReturnedTriples);
                if (numReturnedTriples < 1) {
                    continue;
                }
                const int16_t* matched_triple = NULL;
                bool unique = true;
                /*
                 * basically if we make a query on small angle, all returned catalog triangles have small angle
                 * within toleranceand we want to compare each resulting triangle to affirm the large angle is 
                 * also within tolerance but if there is more than one than we skip this query
                */
                for (const int16_t *l = lowerBoundSearch; l < lowerBoundSearch + numReturnedTriples * 3; l += 3) {
                    float actualLargeAngle = std::max(std::max(Angle(catalog[*(l+1)].spatial-catalog[*l].spatial,
                catalog[*(l+2)].spatial-catalog[*l].spatial), Angle(catalog[*(l+1)].spatial-catalog[*(l+2)].spatial, 
                catalog[*(l+1)].spatial-catalog[*l].spatial)), Angle(catalog[*(l+2)].spatial-catalog[*l].spatial, 
                catalog[*(l+2)].spatial-catalog[*(l+1)].spatial));
                    if (actualLargeAngle > largeAngle - tolerance && actualLargeAngle < largeAngle + tolerance) {
                        if (matched_triple != NULL) {
                            unique = false;
                        }
                        matched_triple = l;
                    }
                }
                // if there was no unique match then find new target triangle
                if (matched_triple == NULL || !unique) {
                    continue;
                }
                // otherwise use random 4th and 5th stars to confirm same triangle
                // TODO
                // ...
                // we have a matching triple
                identified_fast[i] = true;
                identified_fast[j] = true;
                identified_fast[k] = true;
                // TODO figure out which one is which XDXDXDXD (compare angles small medium large)
                StarIdentifier newStar1(i, *matched_triple);
                StarIdentifier newStar2(j, *(matched_triple+1));
                StarIdentifier newStar3(k, *(matched_triple+2));
                // Set identified[i] to value of catalog index of star w most votesr
                identified.push_back(newStar1);
                identified.push_back(newStar2);
                identified.push_back(newStar3);

            }
        }
    }
    return identified;
}

}
