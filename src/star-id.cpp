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
    //make star datastructure that will keep track of the votes?
    for (const Star &i : stars) {
        //convert x and y coordinates to degree differences 
        //ascension or declination = arctan(xpos/(xRes/2/tan(FOV/2)))
        float declination = atan(i.x/(camera.xResolution/2/tan(camera.xFov/2)));
        float rascension = atan(i.y/(camera.xResolution/2/tan(camera.xFov/2)));
        for (const Star &j: stars) {
            if (&i != &j) {
                float secondStarDeclination = atan(j.x/(camera.xResolution/2/tan(camera.xFov/2)));
                float secondStarRascension = atan(i.y/(camera.xResolution/2/tan(camera.xFov/2)));
                float GCD = 2.0*asin(sqrt(pow(sin(abs(declination-secondStarDeclination)/2.0), 2.0)
                         + cos(declination)*cos(secondStarDeclination)*pow(sin(abs(rascension-secondStarRascension)/2.0), 2.0)));
                //give a greter range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = GCD + i.radiusX/2 + j.radiusX/2;
                float upperBoundRange = lowerBoundRange + i.radiusX/2 + j.radiusX/2;
                //if database is a KVectorDatabase
                int *numReturnedPairs; 
                int lowerBoundSearch = database.FindPossiblestarPairsApprox(lowerBoundRange, upperBoundRange, numReturnedPairs);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastruct
                //US voting system 
            }
        }
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
