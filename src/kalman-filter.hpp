#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <array>

#include "attitude-utils.hpp"

namespace lost {

/**
 * Multiplicative extended Kalman filter (MEKF) for attitude and angular rate determination.
 *
 * State: the attitude quaternion `q` (inertial to body) and the rate `w`.
 * The error state is 6-dimensional: a small rotation `theta` applied on the left of `q`
 * (q_true = Exp(theta) * q) and an additive rate error. `P` is its covariance, ordered
 * (theta_x, theta_y, theta_z, w_x, w_y, w_z), in radians and radians per second.
 *
 * `w` is defined by q(t + dt) = Exp(w * dt) * q(t). Since LOST quaternions rotate inertial vectors
 * into the body frame, `w` is the body-frame angular velocity with its sign flipped. This must be
 * kept in mind before feeding in gyro measurements.
 *
 * The only measurement is a full attitude quaternion (e.g. from Davenport Q).
 *
 * The filter does nothing until Reset succeeds: before that, Predict and Update return false and
 * the state is the identity attitude with zero rate and zero covariance. From then on the state
 * is always usable (unit attitude, finite rate and covariance): Reset, Predict and Update return
 * false and change nothing rather than break that.
 */
class AttitudeEKF {
public:
    using Mat6 = std::array<std::array<decimal, 6>, 6>;

    /**
     * @param crossBoresightSigma 1-sigma error across the boresight (body y and z), in radians.
     * @param rollSigma 1-sigma error about the boresight (body x), in radians.
     * @param angularAccelerationSigma Unmodeled angular acceleration density, in rad/s^(3/2).
     * The first two must be positive and the third at least zero, all finite; otherwise Reset
     * always returns false.
     */
    AttitudeEKF(decimal crossBoresightSigma, decimal rollSigma, decimal angularAccelerationSigma);

    /**
     * Start from the newer of two measured attitudes `baseline` seconds apart, and the rate
     * between them: Log(newer * older^-1) / baseline. The initial covariance is what those two
     * measurements imply, so it is only right for a rate obtained this way.
     * @return false if the attitude is not unit, the rate is not finite, the baseline is not
     * positive and finite, the two imply more than half a turn, or the covariance is not
     * representable.
     */
    bool Reset(const Quaternion &attitude, const Vec3 &rate, decimal baseline);

    /**
     * Propagate the state `dt` seconds forward assuming constant rate.
     * @return false before the first successful Reset, if `dt` is negative or not finite, or if
     * the covariance would overflow.
     */
    bool Predict(decimal dt);

    /**
     * Update the state with a measured attitude.
     * @param gate Largest accepted squared Mahalanobis distance between measurement and prediction.
     * If the filter's noise settings match the sensor, that distance is chi-square distributed
     * with 3 degrees of freedom, and a gate of 16.27 rejects 0.1% of good measurements.
     * @return false before the first successful Reset, if the measurement is not unit, or if it
     * fails the gate.
     */
    bool Update(const Quaternion &measured, decimal gate);

    const Quaternion &AttitudeQuaternion() const { return q; }
    const Vec3 &Rate() const { return w; }
    const Mat6 &Covariance() const { return P; }

private:
    Quaternion q;
    Vec3 w;
    Mat6 P;

    /// Measurement noise variances about body x (roll), y and z, in rad^2 (the covariance is diagonal).
    std::array<decimal, 3> R;
    /// Spectral density of the unmodeled angular acceleration, in rad^2/s^3.
    decimal accelerationVariance;
    /// Whether the constructor's arguments can be used at all.
    bool settingsValid;
    /// Whether Reset has succeeded.
    bool initialized;
};

/// How far the squared norm of a quaternion may be from 1 for it to be taken as an attitude.
const decimal kUnitTolerance = DECIMAL(1e-3);

/// Exp: the unit quaternion for rotation vector v, a rotation of |v| radians about v.
Quaternion RotationVectorToQuaternion(const Vec3 &v);
/// Log: the rotation vector of unit quaternion q, taking the rotation of at most pi radians.
/// At exactly pi radians q and -q give opposite vectors, both of them correct.
Vec3 QuaternionToRotationVector(const Quaternion &q);

}  // namespace lost

#endif
