#ifndef ATTITUDE_ESTIMATORS_H
#define ATTITUDE_ESTIMATORS_H

#include <eigen3/Eigen/Eigenvalues>

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-id.hpp"

namespace lost {

/// Necessary matrix information for both the dqm and quest algos
struct DqmQuestHelperMatrices {
    Eigen::Matrix4f K;      /// Davenport matrix
    Eigen::Matrix3f S;      /// B + B^T (where B is the attitude profile matrix)
    Eigen::Vector3f Z;      /// [B23 - B32, B31 - B13, B12 - B21]^T]
    float sigma;            /// tr(B)
};

/**
 * Returns necessary matrix information for calculations in both the dqm and quest algorithms
 */
DqmQuestHelperMatrices DqmQuestHelperMatricesConstructor(const Camera &camera, const Stars &stars, const Catalog &catalog, const StarIdentifiers &starIdentifiers);

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

/**
 * A faster and just as accurate attitude estimator as the Davenport Q algorithm.
 * Historically important as the de facto standard attitude estimation algorithm.
 */
class QuestAlgorithm : public AttitudeEstimationAlgorithm {
public:
    Attitude Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &);
};

}

#endif
