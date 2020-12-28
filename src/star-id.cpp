#include <stdlib.h>

#include "star-id.hpp"

namespace lost {

void DummyStarIdAlgorithm::Go(const unsigned char *database, Stars *stars) const {
    for (Star &star : *stars) {
        if (rand() > RAND_MAX/5) {
            star.identifiedBscIndex = rand() % 9110;
        }
    }
}

void GeometricVotingStarIdAlgorithm::Go(const unsigned char *database, Stars *stars) const {
    // TODO
    ;
}

void PyramidStarIdAlgorithm::Go(const unsigned char *database, Stars *stars) const {
    // TODO
    ;
}

}
