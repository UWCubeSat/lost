#ifndef ATTITUDE_ESTIMATORS_H
#define ATTITUDE_ESTIMATORS_H

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-id.hpp"

namespace lost {

/**
 * An attitude estimation algorithm estimates the orientation of the camera based on identified stars.
 * Subclasses should have a constructor which stores any configuration values in private fields, and then override the Go method to perform the actual attitude estimation.
 */
class AttitudeEstimationAlgorithm {
public:
    /**
     * Actually run the star-id algorithm.
     * Uses the given centroids and star identifiers to come up with an attitude estimate which minimizes error.
     * @todo More detail in return type (eg, whether attitude estimation failed, measure of error)
     */
    virtual Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &) = 0;

    virtual ~AttitudeEstimationAlgorithm() {};
};

/**
 * A slow but reliable attitude estimation algorithm.
 * Requires the Eigen3 library to find eigenvectors of a matrix. On the upside, it finds the optimal attitude estimate, and takes into account information from all identified stars.
 */
class DavenportQAlgorithm : public AttitudeEstimationAlgorithm {
public:
    Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

/**
 * A fast attitude estimator which only takes into account information from two stars.
 * TRIAD is prone to error if either of the selected stars are not good. Even if both selected stars are decent, it's better to use an attitude estimator that takes into account all the stars to better cancel out centroiding noise.
 */
class TriadAlgorithm : public AttitudeEstimationAlgorithm {
public:
    Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

}

#endif
