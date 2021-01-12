#include <stdlib.h>

#include "star-id.hpp"

namespace lost {

StarIdentifiers DummyStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars) const {
    StarIdentifiers result;

    for (int i = 0; i < (int)stars.size(); i++) {
        if (rand() > RAND_MAX/5) {
            result.push_back(StarIdentifier(i, rand() % 2000));
        }
    }

    return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars) const {
    //convert x and y coordinates to degree differences 
    //give a greater range for min-max Query for bigger radius(?) (GreatCircleDistance)
    //us voting system 
    //optimizations? N^2
    //testing, add false stars and see if the accuracy is still good (maybe just 1 or 2 false stars)
    ;
}

StarIdentifiers PyramidStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars) const {
    // TODO
    ;
}

}
