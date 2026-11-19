#include "kalman-filter.hpp"

#include <array>
#include <cmath>

#include "attitude-utils.hpp"
#include "decimal.hpp"

namespace lost {

using Mat6 = AttitudeEKF::Mat6;

Quaternion RotationVectorToQuaternion(const Vec3 &v) {
    decimal angle = v.Magnitude();
    if (angle < DECIMAL(1e-12)) {
        // first order; exact to rounding for angles this small
        return Quaternion(DECIMAL(1.0), v.x / 2, v.y / 2, v.z / 2);
    }
    return Quaternion(v * (DECIMAL(1.0) / angle), angle);
}

Vec3 QuaternionToRotationVector(const Quaternion &q) {
    // q and -q are the same rotation; pick the one with the smaller angle
    decimal sign = q.real < 0 ? DECIMAL(-1.0) : DECIMAL(1.0);
    Vec3 axis = Vec3(q.i, q.j, q.k) * sign;
    decimal sinHalf = axis.Magnitude();
    if (sinHalf < DECIMAL(1e-12)) {
        return axis * DECIMAL(2.0);
    }
    decimal angle = 2 * DECIMAL_ATAN2(sinHalf, q.real * sign);
    return axis * (angle / sinHalf);
}

static Quaternion Normalized(const Quaternion &q) {
    decimal inv = DECIMAL(1.0) / DECIMAL_SQRT(q.real * q.real + q.i * q.i + q.j * q.j + q.k * q.k);
    return Quaternion(q.real * inv, q.i * inv, q.j * inv, q.k * inv);
}

static bool IsFinite(const Vec3 &v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static bool IsFinite(const Mat6 &m) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (!std::isfinite(m[i][j])) return false;
        }
    }
    return true;
}

/**
 * Left Jacobian of the rotation Exp(v): to first order in d,
 * Exp(v + d) = Exp(J * d) * Exp(v).
 */
static void LeftJacobian(const Vec3 &v, decimal jacobian[3][3]) {
    const decimal angle = v.Magnitude();
    // sin(angle)/angle and (1 - cos(angle))/angle, the latter written to avoid cancellation
    decimal sinc = 1, cosc = angle / 2;
    Vec3 axis(0, 0, 0);
    if (angle > DECIMAL(1e-6)) {
        const decimal sinHalf = DECIMAL_SIN(angle / 2);
        sinc = DECIMAL_SIN(angle) / angle;
        cosc = 2 * sinHalf * sinHalf / angle;
        axis = v * (DECIMAL(1.0) / angle);
    }
    const decimal a[3] = {axis.x, axis.y, axis.z};
    const decimal cross[3][3] = {{0, -a[2], a[1]}, {a[2], 0, -a[0]}, {-a[1], a[0], 0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            jacobian[i][j] = (i == j ? sinc : 0) + (1 - sinc) * a[i] * a[j] + cosc * cross[i][j];
        }
    }
}

/// A * B^T
static Mat6 MultiplyTransposed(const Mat6 &a, const Mat6 &b) {
    Mat6 result{};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                result[i][j] += a[i][k] * b[j][k];
            }
        }
    }
    return result;
}

/// A * P * A^T, symmetrized
static Mat6 Congruence(const Mat6 &a, const Mat6 &p) {
    Mat6 ap{};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                ap[i][j] += a[i][k] * p[k][j];
            }
        }
    }
    Mat6 result = MultiplyTransposed(ap, a);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < i; j++) {
            result[i][j] = result[j][i] = (result[i][j] + result[j][i]) / 2;
        }
    }
    return result;
}

AttitudeEKF::AttitudeEKF(decimal crossBoresightSigma, decimal rollSigma, decimal angularAccelerationSigma)
    : q(DECIMAL(1.0), 0, 0, 0), w(0, 0, 0), P{},
      R{rollSigma * rollSigma, crossBoresightSigma * crossBoresightSigma, crossBoresightSigma * crossBoresightSigma},
      accelerationVariance(angularAccelerationSigma * angularAccelerationSigma) { }

bool AttitudeEKF::Reset(const Quaternion &attitude, const Vec3 &rate, decimal baseline) {
    // written so that non-finite values are refused too
    if (!attitude.IsUnit(kUnitTolerance) || !IsFinite(rate) || !(baseline > 0)) {
        return false;
    }
    // rate = (theta2 - theta1) / baseline, attitude = theta2, with independent measurement errors
    Mat6 covariance{};
    for (int i = 0; i < 3; i++) {
        covariance[i][i] = R[i];
        covariance[i][i + 3] = covariance[i + 3][i] = R[i] / baseline;
        covariance[i + 3][i + 3] = 2 * R[i] / (baseline * baseline);
    }
    if (!IsFinite(covariance)) {
        return false;
    }
    q = Normalized(attitude);
    w = rate;
    P = covariance;
    return true;
}

bool AttitudeEKF::Predict(decimal dt) {
    if (!(dt >= 0) || !std::isfinite(dt)) {
        return false;
    }
    Quaternion step = RotationVectorToQuaternion(w * dt);
    Quaternion predicted = Normalized(step * q);

    // Attitude error rotates with the step: Exp(step(theta)) * step = step * Exp(theta).
    // Rate error integrates into attitude error through the same step:
    // Exp((w + dw) * dt) = Exp(J * dw * dt) * step.
    decimal jacobian[3][3];
    LeftJacobian(w * dt, jacobian);
    Mat6 phi{};
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int j = 0; j < 3; j++) {
        Vec3 column = step.Rotate(basis[j]);
        phi[0][j] = column.x;
        phi[1][j] = column.y;
        phi[2][j] = column.z;
        for (int i = 0; i < 3; i++) {
            phi[i][j + 3] = dt * jacobian[i][j];
        }
        phi[j + 3][j + 3] = 1;
    }
    Mat6 covariance = Congruence(phi, P);

    // white angular acceleration, integrated over dt (to first order in the step's angle)
    decimal dt2 = dt * dt;
    for (int i = 0; i < 3; i++) {
        covariance[i][i] += accelerationVariance * dt2 * dt / 3;
        covariance[i][i + 3] += accelerationVariance * dt2 / 2;
        covariance[i + 3][i] += accelerationVariance * dt2 / 2;
        covariance[i + 3][i + 3] += accelerationVariance * dt;
    }

    if (!predicted.IsUnit(kUnitTolerance) || !IsFinite(covariance)) {
        return false;  // overflow; keep the last good state
    }
    q = predicted;
    P = covariance;
    return true;
}

bool AttitudeEKF::Update(const Quaternion &measured, decimal gate) {
    if (!measured.IsUnit(kUnitTolerance)) {
        return false;  // not a rotation (or not finite)
    }
    Vec3 innovationVec = QuaternionToRotationVector(measured * q.Conjugate());
    const decimal innovation[3] = {innovationVec.x, innovationVec.y, innovationVec.z};

    // innovation covariance S = P_theta + R, factored as L * L^T
    decimal s[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            s[i][j] = P[i][j] + (i == j ? R[i] : 0);
        }
    }
    decimal l[3][3] = {};
    for (int j = 0; j < 3; j++) {
        decimal diagonal = s[j][j];
        for (int k = 0; k < j; k++) diagonal -= l[j][k] * l[j][k];
        if (!(diagonal > 0)) return false;
        l[j][j] = DECIMAL_SQRT(diagonal);
        for (int i = j + 1; i < 3; i++) {
            decimal value = s[i][j];
            for (int k = 0; k < j; k++) value -= l[i][k] * l[j][k];
            l[i][j] = value / l[j][j];
        }
    }
    // x <- S^-1 x
    auto solve = [&l](decimal x[3]) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < i; k++) x[i] -= l[i][k] * x[k];
            x[i] /= l[i][i];
        }
        for (int i = 2; i >= 0; i--) {
            for (int k = i + 1; k < 3; k++) x[i] -= l[k][i] * x[k];
            x[i] /= l[i][i];
        }
    };

    decimal weighted[3] = {innovation[0], innovation[1], innovation[2]};
    solve(weighted);
    decimal mahalanobisSq = innovation[0] * weighted[0] + innovation[1] * weighted[1] + innovation[2] * weighted[2];
    if (!(mahalanobisSq <= gate)) return false;

    // gain K = P H^T S^-1, where H = [I 0]; S is symmetric, so row i of K is S^-1 (row i of P_theta)
    decimal gain[6][3];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) gain[i][j] = P[i][j];
        solve(gain[i]);
    }

    decimal correction[6] = {};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) correction[i] += gain[i][j] * innovation[j];
    }
    Quaternion corrected = Normalized(RotationVectorToQuaternion({correction[0], correction[1], correction[2]}) * q);
    Vec3 correctedRate = w + Vec3(correction[3], correction[4], correction[5]);

    // Joseph form: P = (I - K H) P (I - K H)^T + K R K^T
    Mat6 iMinusKh{};
    for (int i = 0; i < 6; i++) {
        iMinusKh[i][i] = 1;
        for (int j = 0; j < 3; j++) iMinusKh[i][j] -= gain[i][j];
    }
    Mat6 covariance = Congruence(iMinusKh, P);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 3; k++) covariance[i][j] += gain[i][k] * R[k] * gain[j][k];
        }
    }

    // That covariance is of the error about the old attitude. Re-express it about the corrected
    // one: Exp(correction + e) = Exp(J * e) * Exp(correction). The rate error is unaffected.
    decimal jacobian[3][3];
    LeftJacobian({correction[0], correction[1], correction[2]}, jacobian);
    Mat6 recenter{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) recenter[i][j] = jacobian[i][j];
        recenter[i + 3][i + 3] = 1;
    }
    covariance = Congruence(recenter, covariance);

    if (!corrected.IsUnit(kUnitTolerance) || !IsFinite(correctedRate) || !IsFinite(covariance)) {
        return false;  // overflow; keep the last good state
    }
    q = corrected;
    w = correctedRate;
    P = covariance;
    return true;
}

}  // namespace lost
