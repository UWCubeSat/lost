#ifndef ATTITUDE_ESTIMATORS_H
#define ATTITUDE_ESTIMATORS_H

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-id.hpp"

namespace lost {

/**
 * @brief
 * @details
 */
class AttitudeEstimationAlgorithm {
public:
    /**
     * @brief
     * @todo More detail in return type (eg, whether attitude estimation failed, measure of error)
     * @return
     */
    virtual Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &) = 0;

    /// @brief
    virtual ~AttitudeEstimationAlgorithm() {};
};

/**
 * @brief
 * @details
 */
class DavenportQAlgorithm : public AttitudeEstimationAlgorithm {
public:
    Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

/**
 * @brief
 * @details
 */
class TriadAlgorithm : public AttitudeEstimationAlgorithm {
public:
    Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

}

#endif
