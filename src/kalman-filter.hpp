#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <array>

#include "attitude-utils.hpp"

namespace lost {

/**
 * Multiplicative extended Kalman filter (MEKF) for attitude and angular rate determination.
 *
 * State: the attitude quaternion `q` (inertial to body, as everywhere in LOST) and the rate `w`.
 * The error state is 6-dimensional: a small rotation `theta` applied on the left of `q`
 * (q_true = Exp(theta) * q) and an additive rate error. `P` is its covariance, ordered
 * (theta_x, theta_y, theta_z, w_x, w_y, w_z), in radians and radians per second.
 *
 * `w` is defined by q(t + dt) = Exp(w * dt) * q(t). Since LOST quaternions rotate inertial vectors
 * into the body frame, `w` is the body-frame angular velocity with its sign flipped. Keep that in
 * mind before feeding in gyro measurements.
 *
 * The only measurement is a full attitude quaternion (e.g. from Davenport Q).
 *
 * The state is always usable: the attitude is a unit quaternion and the rate and covariance are
 * finite. Reset, Predict and Update each return false and leave the state untouched rather than
 * break that, whatever they are given. Until the first successful Reset the state is the identity
 * attitude at zero rate with zero covariance.
 */
class AttitudeEKF {
public:
    using Mat6 = std::array<std::array<decimal, 6>, 6>;

    /**
     * @param crossBoresightSigma 1-sigma error of a measured attitude about the two axes
     * perpendicular to the boresight (body y and z), in radians.
     * @param rollSigma 1-sigma error about the boresight (body x), in radians.
     * @param angularAccelerationSigma Spectral density of the unmodeled angular acceleration, in
     * rad/s^(3/2). Larger values track maneuvers better but give larger search windows.
     */
    AttitudeEKF(decimal crossBoresightSigma, decimal rollSigma, decimal angularAccelerationSigma);

    /**
     * Start the filter from a measured attitude and a rate estimated from two measurements taken
     * `baseline` seconds apart. The initial covariance is what those two measurements imply.
     * @return false if the attitude is not a unit quaternion, the rate is not finite, or the
     * baseline is not positive (or so small that the covariance overflows).
     */
    bool Reset(const Quaternion &attitude, const Vec3 &rate, decimal baseline);

    /**
     * Propagate the state `dt` seconds forward assuming constant rate.
     * @return false if `dt` is negative or not finite, or the covariance would overflow.
     */
    bool Predict(decimal dt);

    /**
     * Update the state with a measured attitude.
     * @param gate Largest accepted squared Mahalanobis distance between measurement and prediction.
     * It is chi-square distributed with 3 degrees of freedom: 16.27 rejects 0.1% of good measurements.
     * @return false if the measurement is not a unit quaternion or fails the gate. Always false if
     * the filter was constructed with a zero or non-finite sigma.
     */
    bool Update(const Quaternion &measured, decimal gate);

    const Quaternion &AttitudeQuaternion() const { return q; }
    const Vec3 &Rate() const { return w; }
    const Mat6 &Covariance() const { return P; }

private:
    Quaternion q;
    Vec3 w;
    Mat6 P;

    std::array<decimal, 3> R;  ///< measurement noise variance, diagonal
    decimal accelerationVariance;
};

/// How far the squared norm of a quaternion may be from 1 for it to be taken as an attitude.
const decimal kUnitTolerance = DECIMAL(1e-3);

/// The rotation Exp(v): angle |v| radians about v.
Quaternion RotationVectorToQuaternion(const Vec3 &v);
/// Inverse of RotationVectorToQuaternion, choosing the rotation of at most pi radians.
Vec3 QuaternionToRotationVector(const Quaternion &q);

}  // namespace lost

#endif
