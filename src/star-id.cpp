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
    // TODO
    ;
}

StarIdentifiers PyramidStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars) const {
    // TODO
    ;
}

}
