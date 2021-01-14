#include <stdlib.h>

#include <math.h>
#include <iostream>
#include "star-id.hpp"
#include "databases.hpp"

namespace lost {

StarIdentifiers DummyStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Camera &camera) const {

    StarIdentifiers result;

    for (int i = 0; i < (int)stars.size(); i++) {
        if (rand() > RAND_MAX/5) {
            result.push_back(StarIdentifier(i, rand() % 2000));
        }
    }

    return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Camera &camera) const {
    KVectorDatabase vectorDatabase(database);
    StarIdentifiers identified;
    //make star datastructure that will keep track of the votes?
    for (int i = 0; i < (int)stars.size(); i++) {
        std::vector<int16_t> votes(10000); // Catalog size
        //convert x and y coordinates to degree differences 
        //ascension or declination = arctan(xpos/(xRes/2/tan(FOV/2)))
        // TODO: xResolution - 1?
        float ra1 = atan((stars[i].x-camera.xResolution/2.0)/(camera.xResolution/2.0/tan(camera.xFov/2)));
        float de1 = atan((stars[i].y-camera.yResolution/2.0)/(camera.xResolution/2.0/tan(camera.xFov/2)));
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                float ra2 = atan((stars[j].x-camera.xResolution/2.0)/(camera.xResolution/2.0/tan(camera.xFov/2)));
                float de2 = atan((stars[j].y-camera.yResolution/2.0)/(camera.xResolution/2.0/tan(camera.xFov/2)));
                float GCD = 2.0*asin(sqrt(pow(sin(abs(de1-de2)/2.0), 2.0)
                         + cos(de1)*cos(de2)*pow(sin(abs(ra1-ra2)/2.0), 2.0)));
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = GCD - tolerance;
                float upperBoundRange = GCD + tolerance;
                //if database is a KVectorDatabase
                long numReturnedPairs; 
                int16_t *lowerBoundSearch = vectorDatabase.FindPossibleStarPairsApprox(lowerBoundRange, upperBoundRange, &numReturnedPairs);
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
        std::cerr << maxVotes << std::endl;
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
    const unsigned char *database, const Stars &stars, const Camera &camera) const {
    // TODO
    ;
}

}
