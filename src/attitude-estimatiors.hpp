#ifndef ATTITUDE_ESTIMATORS_H
#define ATTITUDE_ESTIMATORS_H

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-id.hpp"

namespace lost {

class AttitudeEstimationAlgorithm {
public:
    virtual Attitude Go(const Camera &, const Stars &) = 0;
    virtual ~AttitudeEstimationAlgorithm() { };
};

}

#endif
