#include "kalman-filter.hpp"

#include <array>
#include <cmath>
#include <limits>

#include "attitude-utils.hpp"
#include "decimal.hpp"

namespace lost {

using Mat6 = AttitudeEKF::Mat6;

Quaternion RotationVectorToQuaternion(const Vec3 &v) {
    // Implements Exp(v) = cos(theta/2) + u sin(theta/2), with theta = |v| and u = v / theta
    decimal angle = v.Magnitude();
    // first order; exact to rounding for angles this small
    if (angle < DECIMAL(1e-12)) {
        return Quaternion(DECIMAL(1.0), v.x / 2, v.y / 2, v.z / 2);
    }
    return Quaternion(v * (DECIMAL(1.0) / angle), angle);
}

Vec3 QuaternionToRotationVector(const Quaternion &q) {
    // Implements Log(q) = 2 * vec * arctan(|vec|, w) / |vec|, with w the real part and vec the vector part
    // q and -q are the same rotation; pick the one with the smaller angle
    decimal sign = q.real < 0 ? DECIMAL(-1.0) : DECIMAL(1.0);
    Vec3 axis = Vec3(q.i, q.j, q.k) * sign;
    decimal sinHalf = axis.Magnitude();
    // first order; avoids dividing by zero at the identity
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

/// Below this angle the Jacobian coefficients come from their Taylor series, which is then more
/// accurate than the closed forms (those divide one small number by another).
static const decimal kSeriesAngle = DECIMAL(1e-2);

/// I + a * [v]x + b * [v]x^2, where [v]x is the cross-product matrix of v: [v]x * c = v cross c.
static void SkewPolynomial(const Vec3 &v, decimal a, decimal b, decimal result[3][3]) {
    const decimal c[3] = {v.x, v.y, v.z};
    const decimal cross[3][3] = {{0, -c[2], c[1]}, {c[2], 0, -c[0]}, {-c[1], c[0], 0}};
    const decimal normSq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // [v]x^2 = v * v^T - |v|^2 * I
            result[i][j] = (i == j ? 1 - b * normSq : 0) + a * cross[i][j] + b * c[i] * c[j];
        }
    }
}

/**
 * Left Jacobian of Exp at v: to first order in d, Exp(v + d) = Exp(J * d) * Exp(v).
 *   J = I + (1 - cos(angle)) / angle^2 * [v]x + (angle - sin(angle)) / angle^3 * [v]x^2,
 * with angle = |v|.
 */
static void LeftJacobian(const Vec3 &v, decimal jacobian[3][3]) {
    const decimal angle = v.Magnitude();
    decimal a = DECIMAL(0.5) - angle * angle / 24;
    decimal b = DECIMAL(1.0) / 6 - angle * angle / 120;
    if (angle > kSeriesAngle) {
        const decimal sinHalf = DECIMAL_SIN(angle / 2);
        a = 2 * sinHalf * sinHalf / (angle * angle);  // 1 - cos(angle) = 2 sin^2(angle / 2), without cancellation
        b = (angle - DECIMAL_SIN(angle)) / (angle * angle * angle);
    }
    SkewPolynomial(v, a, b, jacobian);
}

/**
 * Inverse of the left Jacobian at v, for |v| less than a full turn:
 *   J^-1 = I - [v]x / 2 + (1 / angle^2 - cot(angle / 2) / (2 * angle)) * [v]x^2,
 * with angle = |v|.
 */
static void LeftJacobianInverse(const Vec3 &v, decimal inverse[3][3]) {
    const decimal angle = v.Magnitude();
    decimal b = DECIMAL(1.0) / 12 + angle * angle / 720;
    if (angle > kSeriesAngle) {
        b = 1 / (angle * angle) - DECIMAL_COS(angle / 2) / (2 * angle * DECIMAL_SIN(angle / 2));
    }
    SkewPolynomial(v, DECIMAL(-0.5), b, inverse);
}

/**
 * B(t): the integral of the rotation Exp(w * tau) over tau from 0 to t, which is t * J(w * t)
 * with J the left Jacobian. It is how much attitude error a unit of rate error becomes after
 * time t at rate w.
 */
static void IntegratedRotation(const Vec3 &w, decimal t, decimal result[3][3]) {
    LeftJacobian(w * t, result);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) result[i][j] *= t;
    }
}

/// Points and weights of 3-point Gauss-Legendre quadrature on the interval from 0 to 1. The
/// weights are positive, so a sum of positive semidefinite terms formed with them is too.
static const decimal kQuadraturePoints[3] = {DECIMAL(0.1127016653792583), DECIMAL(0.5), DECIMAL(0.8872983346207417)};
static const decimal kQuadratureWeights[3] = {DECIMAL(5.0) / 18, DECIMAL(8.0) / 18, DECIMAL(5.0) / 18};

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
      accelerationVariance(angularAccelerationSigma * angularAccelerationSigma),
      // Written so that NaN is refused too. The squares are checked as well, in case a sigma is
      // so small or so large that its square is not representable.
      settingsValid(crossBoresightSigma > 0 && rollSigma > 0 && angularAccelerationSigma >= 0
                    && R[0] > 0 && std::isfinite(R[0]) && R[1] > 0 && std::isfinite(R[1])
                    && std::isfinite(accelerationVariance)),
      initialized(false) { }

/// The largest turn between the two fixes given to Reset, in radians. A rate taken from two fixes
/// cannot show more than half a turn; the margin allows for rounding.
static const decimal kMaxResetTurn = DECIMAL_M_PI * DECIMAL(1.001);

bool AttitudeEKF::Reset(const Quaternion &attitude, const Vec3 &rate, decimal baseline) {
    // !(x > 0), not x <= 0, so that NaN is refused too
    if (!settingsValid || !attitude.IsUnit(kUnitTolerance) || !IsFinite(rate)
        || !(baseline > 0) || !std::isfinite(baseline)) {
        return false;
    }
    const Vec3 turn = rate * baseline;
    if (!(turn.Magnitude() <= kMaxResetTurn)) {
        return false;
    }

    // The attitude is the newer of two fixes and the rate is Log(newer * older^-1) / baseline.
    // Call the fix errors e1 (newer) and e0 (older): independent, each with covariance R. To first
    // order the attitude error is e1 and the rate error is J^-1 * (e1 - D * e0 + a) / baseline,
    // where D is the rotation of the turn between the fixes, J the left Jacobian at that turn, and
    // a is explained below. So:
    //   attitude covariance R,
    //   attitude-rate covariance R * J^-T / baseline,
    //   rate covariance J^-1 * (R + D * R * D^T + A) * J^-T / baseline^2.
    // With no turn, D and J are the identity and these are R, R / baseline and (2R + A) / baseline^2.
    //
    // a and its covariance A: two fixes give the average rate over the baseline, and it is used as
    // the rate at the newer fix. Angular acceleration between the fixes makes those differ. An
    // impulse of acceleration a time v after the older fix shifts the turn by C(v) times the
    // impulse, with C(v) = B(baseline) - B(baseline - v). For white acceleration of spectral
    // density s, A = s * integral of C * C^T over v from 0 to baseline, by Gauss-Legendre
    // quadrature. With no turn A = s * baseline^3 / 3, so the rate variance gains s * baseline / 3.
    decimal scaled[3][3];  // J^-1 / baseline
    LeftJacobianInverse(turn, scaled);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) scaled[i][j] /= baseline;
    }
    const Quaternion turnRotation = RotationVectorToQuaternion(turn);
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    decimal d[3][3];
    for (int j = 0; j < 3; j++) {
        const Vec3 column = turnRotation.Rotate(basis[j]);
        d[0][j] = column.x;
        d[1][j] = column.y;
        d[2][j] = column.z;
    }
    decimal both[3][3];  // R + D * R * D^T + A
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            both[i][j] = i == j ? R[i] : 0;
            for (int k = 0; k < 3; k++) both[i][j] += d[i][k] * R[k] * d[j][k];
        }
    }
    decimal whole[3][3];  // B(baseline)
    IntegratedRotation(rate, baseline, whole);
    for (int point = 0; point < 3; point++) {
        decimal c[3][3];  // C(v) at this quadrature point
        IntegratedRotation(rate, baseline * (1 - kQuadraturePoints[point]), c);
        const decimal weight = accelerationVariance * baseline * kQuadratureWeights[point];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) c[i][j] = whole[i][j] - c[i][j];
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) both[i][j] += weight * (c[i][k] * c[j][k]);
            }
        }
    }

    Mat6 covariance{};
    for (int i = 0; i < 3; i++) {
        covariance[i][i] = R[i];
        for (int j = 0; j < 3; j++) {
            covariance[i][j + 3] = covariance[j + 3][i] = R[i] * scaled[j][i];
            decimal sum = 0;
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) sum += scaled[i][k] * both[k][l] * scaled[j][l];
            }
            covariance[i + 3][j + 3] = sum;
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < i; j++) {
            const decimal mean = (covariance[i + 3][j + 3] + covariance[j + 3][i + 3]) / 2;
            covariance[i + 3][j + 3] = covariance[j + 3][i + 3] = mean;
        }
        // A baseline so long that the rate variance underflows would leave a covariance that is
        // not positive definite.
        if (!(covariance[i + 3][i + 3] >= std::numeric_limits<decimal>::min())) {
            return false;
        }
    }
    // a baseline so short that the covariance overflows
    if (!IsFinite(covariance)) {
        return false;
    }
    q = Normalized(attitude);
    w = rate;
    P = covariance;
    initialized = true;
    return true;
}

bool AttitudeEKF::Predict(decimal dt) {
    // refuse a filter that was never started, and a negative or non-finite step
    if (!initialized || !(dt >= 0) || !std::isfinite(dt)) {
        return false;
    }
    // State: q <- Exp(w * dt) * q. The rate stays as it is.
    Quaternion step = RotationVectorToQuaternion(w * dt);
    Quaternion predicted = Normalized(step * q);

    // Transition matrix phi of the error, ordered attitude then rate:
    //   attitude block: the step's rotation matrix. The attitude error rotates with the step,
    //     Exp(step(theta)) * step = step * Exp(theta), exactly.
    //   rate-to-attitude block: dt * J, with J the left Jacobian at w * dt. A rate error
    //     integrates into attitude error: to first order in dw,
    //     Exp((w + dw) * dt) = Exp(J * dw * dt) * step.
    //   rate block: the identity.
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
    // P <- phi * P * phi^T
    Mat6 covariance = Congruence(phi, P);

    // Process noise of white angular acceleration with spectral density s = accelerationVariance.
    // An impulse of acceleration with time tau left in the step adds B(tau) times the impulse to
    // the attitude error and the impulse itself to the rate error, so the noise covariance is
    //   s * integral over tau from 0 to dt of [B * B^T, B; B^T, I].
    // By Gauss-Legendre quadrature: each point adds a positive multiple of [B; I] * [B; I]^T, so
    // the sum is positive semidefinite however far the step turns. It is exact with no turn
    // (s * dt^3 / 3, s * dt^2 / 2, s * dt) and within 0.1% of each block's scale for a turn of up
    // to half a radian; the tests give the error for larger turns.
    for (int point = 0; point < 3; point++) {
        decimal b[3][3];
        IntegratedRotation(w, dt * kQuadraturePoints[point], b);
        const decimal weight = accelerationVariance * dt * kQuadratureWeights[point];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) covariance[i][j] += weight * (b[i][k] * b[j][k]);
                covariance[i][j + 3] += weight * b[i][j];
                covariance[j + 3][i] += weight * b[i][j];
            }
            covariance[i + 3][i + 3] += weight;
        }
    }

    // overflow; keep the last good state
    if (!predicted.IsUnit(kUnitTolerance) || !IsFinite(covariance)) {
        return false;
    }
    q = predicted;
    P = covariance;
    return true;
}

bool AttitudeEKF::Update(const Quaternion &measured, decimal gate) {
    // refuse a filter that was never started, and a measurement that is not a rotation (or not finite)
    if (!initialized || !measured.IsUnit(kUnitTolerance)) {
        return false;
    }
    // Innovation: Log(measured * q^-1), the rotation from the predicted attitude to the measured one.
    Vec3 innovationVec = QuaternionToRotationVector(measured * q.Conjugate());
    const decimal innovation[3] = {innovationVec.x, innovationVec.y, innovationVec.z};

    // Innovation covariance S = P_theta + R, where P_theta is the attitude block of P.
    decimal s[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            s[i][j] = P[i][j] + (i == j ? R[i] : 0);
        }
    }
    // Cholesky factorization S = L * L^T. Fails if S is not positive definite.
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
    // x <- S^-1 * x, by forward then back substitution with L
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

    // Gate on the squared Mahalanobis distance innovation^T * S^-1 * innovation.
    decimal weighted[3] = {innovation[0], innovation[1], innovation[2]};
    solve(weighted);
    decimal mahalanobisSq = innovation[0] * weighted[0] + innovation[1] * weighted[1] + innovation[2] * weighted[2];
    if (!(mahalanobisSq <= gate)) return false;

    // Gain K = P * H^T * S^-1, where H = [I 0] picks out the attitude error.
    // S is symmetric, so row i of K is S^-1 * (the first three entries of row i of P).
    decimal gain[6][3];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) gain[i][j] = P[i][j];
        solve(gain[i]);
    }

    // Correction K * innovation, ordered attitude then rate.
    decimal correction[6] = {};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) correction[i] += gain[i][j] * innovation[j];
    }
    // Apply it: q <- Exp(attitude correction) * q, w <- w + rate correction.
    Quaternion corrected = Normalized(RotationVectorToQuaternion({correction[0], correction[1], correction[2]}) * q);
    Vec3 correctedRate = w + Vec3(correction[3], correction[4], correction[5]);

    // Joseph form: P <- (I - K H) * P * (I - K H)^T + K * R * K^T
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
    // one: P <- G * P * G^T with G = diag(J, I), J the left Jacobian at the attitude correction,
    // since to first order in e, Exp(correction + e) = Exp(J * e) * Exp(correction). The rate
    // error is unaffected.
    decimal jacobian[3][3];
    LeftJacobian({correction[0], correction[1], correction[2]}, jacobian);
    Mat6 recenter{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) recenter[i][j] = jacobian[i][j];
        recenter[i + 3][i + 3] = 1;
    }
    covariance = Congruence(recenter, covariance);

    // overflow; keep the last good state
    if (!corrected.IsUnit(kUnitTolerance) || !IsFinite(correctedRate) || !IsFinite(covariance)) {
        return false;
    }
    q = corrected;
    w = correctedRate;
    P = covariance;
    return true;
}

}  // namespace lost
