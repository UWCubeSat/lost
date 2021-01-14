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
        std::vector<int16_t> votes(catalog.size());
        float ra1, de1;
        camera.CoordinateAngles({ stars[i].x, stars[i].y }, &ra1, &de1);
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                float ra2, de2;
                camera.CoordinateAngles({ stars[j].x, stars[j].y }, &ra2, &de2);
                float gcd = GreatCircleDistance(ra1, de1, ra2, de2);
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = gcd - tolerance;
                float upperBoundRange = gcd + tolerance;
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
        for (int v = 0; v < (int)votes.size(); v++) {
            if (votes[v] > maxVotes) {
                maxVotes = votes[v];
                indexOfMax = v;
            }
        }
        //starIndex = i, catalog index = indexOfMax
        StarIdentifier newStar(i, indexOfMax);
        // Set identified[i] to value of catalog index of star w most votesr
        identified.push_back(newStar);
    }
    //optimizations? N^2
    //testing, add false stars and see if the accuracy is still good (maybe just 1 or 2 false stars)
    return identified;
}

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {
    // TODO
    ;
}

}
