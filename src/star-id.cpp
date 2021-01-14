#include <stdlib.h>

#include <math.h>
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
    for (auto i = 0; i < stars.size(); i++) {
        std::vector<int16_t> votes(10000); // Catalog size
        //convert x and y coordinates to degree differences 
        //ascension or declination = arctan(xpos/(xRes/2/tan(FOV/2)))
        float declination = atan(stars[i].x/(camera.xResolution/2/tan(camera.xFov/2)));
        float rascension = atan(stars[i].y/(camera.xResolution/2/tan(camera.xFov/2)));
        for (const Star &j: stars) {
            if (&stars[i] != &j) {
                float secondStarDeclination = atan(j.x/(camera.xResolution/2/tan(camera.xFov/2)));
                float secondStarRascension = atan(stars[i].y/(camera.xResolution/2/tan(camera.xFov/2)));
                float GCD = 2.0*asin(sqrt(pow(sin(abs(declination-secondStarDeclination)/2.0), 2.0)
                         + cos(declination)*cos(secondStarDeclination)*pow(sin(abs(rascension-secondStarRascension)/2.0), 2.0)));
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = GCD - tolerance;
                float upperBoundRange = GCD + tolerance;
                //if database is a KVectorDatabase
                int numReturnedPairs; 
                int16_t *lowerBoundSearch = vectorDatabase.FindPossibleStarPairsApprox(lowerBoundRange, upperBoundRange, &numReturnedPairs);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastructure
                for (auto k = 0; k < numReturnedPairs; k++) {
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
        for (auto v = 0; v < votes.size(); v++) {
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
}

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Camera &camera) const {
    // TODO
    ;
}

}
