#include <stdlib.h>

#include <math.h>
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
    KVectorDatabase vectorDatabase(database);
    StarIdentifiers identified;
    for (int i = 0; i < (int)stars.size(); i++) {  
        std::vector<int16_t> votes(catalog.size(), 0);
        // TODO: store spatial coordinates in catalog to avoid this conversion
        Vec3 iSpatial = camera.CameraToSpatial({ stars[i].x, stars[i].y });
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                Vec3 jSpatial = camera.CameraToSpatial({ stars[j].x, stars[j].y });
                float greatCircleDistance = AngleUnit(iSpatial, jSpatial);
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = greatCircleDistance - tolerance;
                float upperBoundRange = greatCircleDistance + tolerance;
                //if database is a KVectorDatabase
                long numReturnedPairs; 
                int16_t *lowerBoundSearch = vectorDatabase.FindPossibleStarPairsApprox(
                    lowerBoundRange, upperBoundRange, &numReturnedPairs);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastructure
                for (long k = 0; k < numReturnedPairs; k++) {
                    int16_t first = *(lowerBoundSearch + 2 * k);
                    int16_t second = *(lowerBoundSearch + 2 * k + 1);
                    votes[first]++;
                    votes[second]++;
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
        printf("Max votes: %d\n", maxVotes);
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
    int16_t maxVotes = verificationVotes[0];
    for (int v = 1; v < (int)verificationVotes.size(); v++) {
        if (verificationVotes[v] > maxVotes) {
            maxVotes = verificationVotes[v];
        }
    }

    // If the stars are within a certain range of the maximal number of votes, 
    // we consider it correct.
    // maximal votes = maxVotes
    StarIdentifiers verified;
    int16_t thresholdVotes = 0.75 * maxVotes;
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

}
