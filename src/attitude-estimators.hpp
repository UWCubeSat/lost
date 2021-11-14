#ifndef ATTITUDE_ESTIMATORS_H
#define ATTITUDE_ESTIMATORS_H

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-id.hpp"

namespace lost {

class AttitudeEstimationAlgorithm {
public:
    // TODO: more detail in return type (eg, whether attitude estimation failed, measure of error)
    virtual Quaternion Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &) = 0;
    virtual ~AttitudeEstimationAlgorithm() { };
};

class DavenportQAlgorithm : public AttitudeEstimationAlgorithm {
    public:
        Quaternion Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

}

#endif
