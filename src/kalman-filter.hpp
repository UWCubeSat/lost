#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <array>

#include "attitude-utils.hpp"

namespace lost {

/**
 * Multiplicative extended Kalman filter (MEKF) for attitude and angular rate.
 *
 * State: the attitude quaternion `q` (inertial to body) and the rate `w`.
 * Error state: a small rotation `theta` on the left of `q` (q_true = Exp(theta) * q) and an
 * additive rate error. `P` is its covariance, ordered (theta_x, theta_y, theta_z, w_x, w_y, w_z),
 * in radians and radians per second.
 *
 * `w` is defined by q(t + dt) = Exp(w * dt) * q(t). LOST quaternions rotate inertial vectors into
 * the body frame, so `w` is the body-frame angular velocity with its sign flipped. This must be
 * kept in mind before feeding in gyro measurements.
 *
 * The only measurement is a full attitude quaternion (e.g. from Davenport Q).
 *
 * Nothing happens until Reset succeeds. After that the state is always usable: unit attitude,
 * finite rate, positive definite covariance. Reset, Predict and Update return false and change
 * nothing rather than break that. If one keeps returning false, start again with Reset.
 *
 * Float mode cannot hold attitude and rate errors that are almost perfectly correlated. That
 * happens after about a hundred Reset baselines of Predict with no Update and almost no process
 * noise. Predict then returns false.
 */
class AttitudeEKF {
public:
    using Mat6 = std::array<std::array<decimal, 6>, 6>;

    /**
     * @param crossBoresightSigma 1-sigma error across the boresight (body y and z), in radians.
     * @param rollSigma 1-sigma error about the boresight (body x), in radians.
     * @param angularAccelerationSigma Unmodeled angular acceleration density, in rad/s^(3/2).
     * All finite, the first two positive, the third not negative. Otherwise Reset returns false.
     */
    AttitudeEKF(decimal crossBoresightSigma, decimal rollSigma, decimal angularAccelerationSigma);

    /**
     * Start from the newer of two fixes `baseline` seconds apart and the rate between them,
     * Log(newer * older^-1) / baseline. The covariance assumes the rate was found this way.
     * @return false if the attitude is not unit, the rate is not finite, the baseline is not
     * positive and finite, the fixes are more than a quarter turn apart, or the covariance
     * cannot be held.
     */
    bool Reset(const Quaternion &attitude, const Vec3 &rate, decimal baseline);

    /**
     * Propagate the state `dt` seconds at constant rate. One long call and many short calls
     * give nearly the same covariance (tested: one step against 1000, turns up to 20 rad, within
     * 0.1% of sqrt(P_ii * P_jj)). Very many tiny steps add rounding error, most in float mode.
     * @return false before Reset, if `dt` is negative or not finite, if the step turns more than
     * 64 rad, or if the covariance overflows or cannot be held.
     */
    bool Predict(decimal dt);

    /**
     * Update the state with a measured attitude.
     * @param gate Largest accepted squared Mahalanobis distance. With correct noise settings it
     * is chi-square with 3 degrees of freedom: 16.27 rejects 0.1% of good measurements.
     * @return false before Reset, if the measurement is not unit, if it fails the gate, or if
     * the covariance cannot be held.
     */
    bool Update(const Quaternion &measured, decimal gate);

    const Quaternion &AttitudeQuaternion() const { return q; }
    const Vec3 &Rate() const { return w; }
    const Mat6 &Covariance() const { return P; }

private:
    Quaternion q;
    Vec3 w;
    Mat6 P;

    /// Measurement noise variances about body x (roll), y and z, in rad^2.
    std::array<decimal, 3> R;
    /// Angular acceleration spectral density, in rad^2/s^3.
    decimal accelerationVariance;
    /// True if the constructor's sigmas are usable.
    bool settingsValid;
    /// True once Reset has succeeded.
    bool initialized;
};

/// Largest | |q|^2 - 1 | for a quaternion to count as an attitude.
const decimal kUnitTolerance = DECIMAL(1e-3);

/// Exp: unit quaternion for rotation vector v, |v| radians about v.
Quaternion RotationVectorToQuaternion(const Vec3 &v);
/// Log: rotation vector of unit quaternion q, at most pi radians.
/// At exactly pi, q and -q give opposite vectors. Both are correct.
Vec3 QuaternionToRotationVector(const Quaternion &q);

}  // namespace lost

#endif
