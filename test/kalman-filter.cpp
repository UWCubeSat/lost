#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <utility>

#include <catch.hpp>

#include "attitude-utils.hpp"
#include "kalman-filter.hpp"

using namespace lost; // NOLINT

namespace {

const decimal kArcsec = DECIMAL(M_PI) / (180 * 3600);

decimal AngleBetween(const Quaternion &a, const Quaternion &b) {
    return QuaternionToRotationVector(a * b.Conjugate()).Magnitude();
}

decimal NormSq(const Quaternion &q) {
    return q.real * q.real + q.i * q.i + q.j * q.j + q.k * q.k;
}

/// Standard normal samples from a generator whose output is fixed by the C++ standard.
class Gaussian {
public:
    explicit Gaussian(unsigned seed) : generator(seed) { }
    decimal operator()() {
        // Box-Muller; (0, 1] uniforms
        const double u1 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
        const double u2 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
        return static_cast<decimal>(std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2));
    }
private:
    std::mt19937 generator;
};

/// e^T * P^-1 * e in double precision, by Gaussian elimination with partial pivoting.
double Nees(const AttitudeEKF::Mat6 &covariance, const double error[6]) {
    double a[6][7];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) a[i][j] = static_cast<double>(covariance[i][j]);
        a[i][6] = error[i];
    }
    for (int col = 0; col < 6; col++) {
        int pivot = col;
        for (int row = col + 1; row < 6; row++) {
            if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;
        }
        for (int j = 0; j < 7; j++) std::swap(a[col][j], a[pivot][j]);
        for (int row = col + 1; row < 6; row++) {
            const double factor = a[row][col] / a[col][col];
            for (int j = col; j < 7; j++) a[row][j] -= factor * a[col][j];
        }
    }
    double solution[6];
    for (int i = 5; i >= 0; i--) {
        double value = a[i][6];
        for (int j = i + 1; j < 6; j++) value -= a[i][j] * solution[j];
        solution[i] = value / a[i][i];
    }
    double result = 0;
    for (int i = 0; i < 6; i++) result += error[i] * solution[i];
    return result;
}

/**
 * Whether a covariance is positive definite. Judged on its correlation matrix, so that the
 * different scales of attitude and rate do not matter: Cholesky in double precision, with every
 * pivot at least `tolerance`. (A pivot is the share of a component's variance that the
 * components before it do not explain.)
 */
bool IsPositiveDefinite(const AttitudeEKF::Mat6 &covariance, double tolerance) {
    double l[6][6] = {};
    for (int j = 0; j < 6; j++) {
        if (!(covariance[j][j] > 0)) return false;
        for (int i = j; i < 6; i++) {
            if (!(covariance[i][i] > 0)) return false;
            double value = static_cast<double>(covariance[i][j])
                / std::sqrt(static_cast<double>(covariance[i][i]) * static_cast<double>(covariance[j][j]));
            for (int k = 0; k < j; k++) value -= l[i][k] * l[j][k];
            if (i == j) {
                if (!(value >= tolerance)) return false;
                l[j][j] = std::sqrt(value);
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    return true;
}

/// A spacecraft whose rate wanders with white angular acceleration, and a star tracker on it.
class NoisyTruth {
public:
    NoisyTruth(const Quaternion &attitude, const Vec3 &rate, decimal crossSigma, decimal rollSigma,
               decimal accelerationSigma, Gaussian *gaussian)
        : attitude(attitude), rate(rate), crossSigma(crossSigma), rollSigma(rollSigma),
          accelerationSigma(accelerationSigma), gaussian(gaussian) { }

    /// Move on by dt: the rate change over dt and the attitude change it causes.
    void Advance(decimal dt) {
        decimal rateStep[3], angleStep[3];
        for (int i = 0; i < 3; i++) {
            const decimal n1 = (*gaussian)(), n2 = (*gaussian)();
            rateStep[i] = accelerationSigma * DECIMAL_SQRT(dt) * n1;
            angleStep[i] = accelerationSigma * dt * DECIMAL_SQRT(dt) * (n1 / 2 + n2 / DECIMAL_SQRT(DECIMAL(12.0)));
        }
        attitude = RotationVectorToQuaternion(rate * dt + Vec3(angleStep[0], angleStep[1], angleStep[2])) * attitude;
        rate = rate + Vec3(rateStep[0], rateStep[1], rateStep[2]);
    }

    Quaternion Measure() {
        const Vec3 noise(rollSigma * (*gaussian)(), crossSigma * (*gaussian)(), crossSigma * (*gaussian)());
        return RotationVectorToQuaternion(noise) * attitude;
    }

    Quaternion attitude;
    Vec3 rate;

private:
    decimal crossSigma, rollSigma, accelerationSigma;
    Gaussian *gaussian;
};

/// The rotation matrix of q, in double precision.
void RotationMatrix(const Quaternion &q, double matrix[3][3]) {
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int j = 0; j < 3; j++) {
        const Vec3 column = q.Rotate(basis[j]);
        matrix[0][j] = static_cast<double>(column.x);
        matrix[1][j] = static_cast<double>(column.y);
        matrix[2][j] = static_cast<double>(column.z);
    }
}

}  // namespace

TEST_CASE("Rotation vector round trip, and q and -q give the same rotation vector", "[kalman]") {
    const Vec3 v(DECIMAL(0.3), DECIMAL(-0.2), DECIMAL(0.1));
    const Quaternion q = RotationVectorToQuaternion(v);
    const Vec3 back = QuaternionToRotationVector(q);
    CHECK(back.x == Approx(v.x).margin(1e-6));
    CHECK(back.y == Approx(v.y).margin(1e-6));
    CHECK(back.z == Approx(v.z).margin(1e-6));
    const Vec3 negated = QuaternionToRotationVector(Quaternion(-q.real, -q.i, -q.j, -q.k));
    CHECK(negated.x == Approx(v.x).margin(1e-6));
}

TEST_CASE("EKF prediction with a constant rate is the exact rotation", "[kalman]") {
    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(0.01));
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.2), DECIMAL(0.3));
    const Vec3 rate(DECIMAL(0.02), DECIMAL(-0.05), DECIMAL(0.01));
    ekf.Reset(start, rate, DECIMAL(0.1));
    for (int i = 0; i < 10; i++) {
        ekf.Predict(DECIMAL(0.1));
    }
    // ten steps compose to one rotation of the same rate (tolerance allows for float mode)
    const Quaternion expected = RotationVectorToQuaternion(rate * DECIMAL(1.0)) * start;
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), expected) < kArcsec / 2);
}

TEST_CASE("EKF tracks a constant rotation from noiseless measurements", "[kalman]") {
    const decimal dt = DECIMAL(0.1);
    const Vec3 trueRate(DECIMAL(0.001), DECIMAL(0.03), DECIMAL(-0.01));
    const Quaternion start = SphericalToQuaternion(DECIMAL(2.0), DECIMAL(-0.4), DECIMAL(1.0));
    auto truth = [&](int frame) { return RotationVectorToQuaternion(trueRate * (frame * dt)) * start; };

    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(0.01));
    // rate deliberately a bit off, as if learned from noisy fixes
    ekf.Reset(truth(1), trueRate * DECIMAL(1.05), dt);
    for (int frame = 2; frame < 50; frame++) {
        ekf.Predict(dt);
        REQUIRE(ekf.Update(truth(frame), DECIMAL(16.27)));
    }
    ekf.Predict(dt);
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), truth(50)) < kArcsec / 10);
    CHECK((ekf.Rate() - trueRate).Magnitude() < DECIMAL(1e-5));
}

TEST_CASE("EKF rejects a measurement that contradicts the prediction, leaving its state alone", "[kalman]") {
    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(0.01));
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.0));
    ekf.Reset(start, Vec3(0, 0, 0), DECIMAL(0.1));
    ekf.Predict(DECIMAL(0.1));
    const Quaternion before = ekf.AttitudeQuaternion();

    const Quaternion farOff = RotationVectorToQuaternion(Vec3(0, 600 * kArcsec, 0)) * start;
    CHECK_FALSE(ekf.Update(farOff, DECIMAL(16.27)));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), before) == 0);

    // not attitudes at all
    CHECK_FALSE(ekf.Update(Quaternion(0, 0, 0, 0), DECIMAL(16.27)));
    CHECK_FALSE(ekf.Update(Quaternion(std::numeric_limits<decimal>::quiet_NaN(), 0, 0, 0), DECIMAL(16.27)));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), before) == 0);

    // q and -q are the same attitude
    const Quaternion close = RotationVectorToQuaternion(Vec3(0, 2 * kArcsec, 0)) * start;
    CHECK(ekf.Update(Quaternion(-close.real, -close.i, -close.j, -close.k), DECIMAL(16.27)));
}

TEST_CASE("EKF refuses inputs that would corrupt its state", "[kalman]") {
    const decimal nan = std::numeric_limits<decimal>::quiet_NaN();
    const decimal inf = std::numeric_limits<decimal>::infinity();
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.0));
    const Vec3 rate(DECIMAL(0.01), 0, 0);

    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(0.01));
    REQUIRE(ekf.Reset(start, rate, DECIMAL(0.1)));
    REQUIRE(ekf.Predict(DECIMAL(0.1)));
    const Quaternion before = ekf.AttitudeQuaternion();
    const AttitudeEKF::Mat6 covarianceBefore = ekf.Covariance();

    CHECK_FALSE(ekf.Reset(Quaternion(0, 0, 0, 0), rate, DECIMAL(0.1)));
    CHECK_FALSE(ekf.Reset(Quaternion(2, 0, 0, 0), rate, DECIMAL(0.1)));
    CHECK_FALSE(ekf.Reset(Quaternion(nan, 0, 0, 0), rate, DECIMAL(0.1)));
    CHECK_FALSE(ekf.Reset(start, Vec3(nan, 0, 0), DECIMAL(0.1)));
    CHECK_FALSE(ekf.Reset(start, Vec3(0, inf, 0), DECIMAL(0.1)));
    CHECK_FALSE(ekf.Reset(start, rate, 0));
    CHECK_FALSE(ekf.Reset(start, rate, DECIMAL(-0.1)));
    CHECK_FALSE(ekf.Reset(start, rate, nan));
    CHECK_FALSE(ekf.Reset(start, rate, inf));
    CHECK_FALSE(ekf.Reset(start, rate, std::numeric_limits<decimal>::min()));  // covariance overflows
    // so long a baseline that the rate variance underflows, which would leave the covariance indefinite
    CHECK_FALSE(ekf.Reset(start, Vec3(0, 0, 0), DECIMAL_SQRT(std::numeric_limits<decimal>::max())));
    // more than a quarter turn between the two fixes
    CHECK_FALSE(ekf.Reset(start, Vec3(DECIMAL(4.0), 0, 0), DECIMAL(1.0)));
    CHECK_FALSE(ekf.Reset(start, Vec3(DECIMAL(1.6), 0, 0), DECIMAL(1.0)));
    CHECK(ekf.Reset(start, Vec3(DECIMAL(1.5), 0, 0), DECIMAL(1.0)));
    REQUIRE(ekf.Reset(start, rate, DECIMAL(0.1)));
    REQUIRE(ekf.Predict(DECIMAL(0.1)));

    CHECK_FALSE(ekf.Predict(DECIMAL(-0.1)));
    CHECK_FALSE(ekf.Predict(nan));
    CHECK_FALSE(ekf.Predict(inf));
    CHECK_FALSE(ekf.Predict(std::numeric_limits<decimal>::max() / 2));  // covariance overflows

    CHECK_FALSE(ekf.Update(before, nan));

    // none of that changed anything
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), before) == 0);
    CHECK((ekf.Rate() - rate).Magnitude() == 0);
    CHECK(ekf.Covariance() == covarianceBefore);

    // and the filter still works
    CHECK(ekf.Predict(0));
    CHECK(ekf.Update(before, DECIMAL(16.27)));
}

TEST_CASE("EKF with unusable noise settings never starts", "[kalman]") {
    const decimal nan = std::numeric_limits<decimal>::quiet_NaN();
    const decimal inf = std::numeric_limits<decimal>::infinity();
    const decimal tiny = std::numeric_limits<decimal>::min();  // its square underflows to zero
    const decimal cross = 5 * kArcsec, roll = 30 * kArcsec, acceleration = DECIMAL(0.01);
    const Quaternion identity(1, 0, 0, 0);
    const decimal settings[][3] = {
        {0, roll, acceleration}, {cross, 0, acceleration}, {-cross, roll, acceleration}, {cross, -roll, acceleration},
        {nan, roll, acceleration}, {cross, nan, acceleration}, {inf, roll, acceleration}, {cross, inf, acceleration},
        {tiny, roll, acceleration}, {cross, tiny, acceleration},
        {cross, roll, -acceleration}, {cross, roll, nan}, {cross, roll, inf},
    };
    for (const decimal *setting : settings) {
        INFO("sigmas " << setting[0] << ", " << setting[1] << ", " << setting[2]);
        AttitudeEKF ekf(setting[0], setting[1], setting[2]);
        CHECK_FALSE(ekf.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
        CHECK_FALSE(ekf.Predict(DECIMAL(0.1)));
        CHECK_FALSE(ekf.Update(identity, DECIMAL(16.27)));
        CHECK(ekf.AttitudeQuaternion().IsUnit(kUnitTolerance));
    }

    // no unmodeled acceleration at all is a usable setting
    AttitudeEKF steady(cross, roll, 0);
    CHECK(steady.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
    CHECK(steady.Predict(DECIMAL(0.1)));
    CHECK(steady.Update(identity, DECIMAL(16.27)));
}

TEST_CASE("EKF does nothing until Reset succeeds", "[kalman]") {
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.0));
    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(0.01));
    CHECK_FALSE(ekf.Predict(DECIMAL(0.1)));
    CHECK_FALSE(ekf.Update(Quaternion(1, 0, 0, 0), DECIMAL(16.27)));
    // a refused Reset does not start it either
    CHECK_FALSE(ekf.Reset(Quaternion(2, 0, 0, 0), Vec3(0, 0, 0), DECIMAL(0.1)));
    CHECK_FALSE(ekf.Predict(DECIMAL(0.1)));
    CHECK_FALSE(ekf.Update(Quaternion(1, 0, 0, 0), DECIMAL(16.27)));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), Quaternion(1, 0, 0, 0)) == 0);
    CHECK(ekf.Rate().Magnitude() == 0);
    CHECK(ekf.Covariance() == AttitudeEKF::Mat6{});

    CHECK(ekf.Reset(start, Vec3(DECIMAL(0.01), 0, 0), DECIMAL(0.1)));
    CHECK(ekf.Predict(DECIMAL(0.1)));
    CHECK(ekf.Update(ekf.AttitudeQuaternion(), DECIMAL(16.27)));
}

TEST_CASE("EKF covariance prediction matches finite differences of the state prediction", "[kalman]") {
    // a large step, so that any shortcut that only holds for small rotations shows
    const Quaternion start = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));
    const Vec3 rate(DECIMAL(0.3), DECIMAL(-0.2), DECIMAL(1.5));
    const decimal dt = DECIMAL(1.0);
    const decimal baseline = DECIMAL(0.5);
    const decimal accelerationSigma = DECIMAL(1e-5);
    const decimal epsilon = DECIMAL(1e-2);

    auto predicted = [&](const Quaternion &attitude, const Vec3 &startRate) {
        AttitudeEKF ekf(DECIMAL(1e-3), DECIMAL(3e-3), accelerationSigma);
        REQUIRE(ekf.Reset(attitude, startRate, baseline));
        REQUIRE(ekf.Predict(dt));
        return ekf;
    };
    const AttitudeEKF nominal = predicted(start, rate);
    AttitudeEKF atReset(DECIMAL(1e-3), DECIMAL(3e-3), accelerationSigma);
    REQUIRE(atReset.Reset(start, rate, baseline));
    const AttitudeEKF::Mat6 before = atReset.Covariance();

    // state transition by central differences: column j is the change of the predicted error state
    // per unit of error component j before the prediction
    AttitudeEKF::Mat6 phi{};
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int j = 0; j < 6; j++) {
        Quaternion attitudes[2] = {start, start};
        Vec3 rates[2] = {rate, rate};
        for (int side = 0; side < 2; side++) {
            const decimal step = side == 0 ? epsilon : -epsilon;
            if (j < 3) {
                attitudes[side] = RotationVectorToQuaternion(basis[j] * step) * start;
            } else {
                rates[side] = rate + basis[j - 3] * step;
            }
        }
        const AttitudeEKF plus = predicted(attitudes[0], rates[0]);
        const AttitudeEKF minus = predicted(attitudes[1], rates[1]);
        const Vec3 attitudeChange = QuaternionToRotationVector(
            plus.AttitudeQuaternion() * minus.AttitudeQuaternion().Conjugate()) * (1 / (2 * epsilon));
        const Vec3 rateChange = (plus.Rate() - minus.Rate()) * (1 / (2 * epsilon));
        const decimal column[6] = {attitudeChange.x, attitudeChange.y, attitudeChange.z,
                                   rateChange.x, rateChange.y, rateChange.z};
        for (int i = 0; i < 6; i++) phi[i][j] = column[i];
    }

    // a rate error perpendicular to the turn axis does not simply add dt times itself to the attitude
    CHECK(std::abs(phi[0][3] - dt) > DECIMAL(0.05));

    // expected covariance: phi * before * phi^T (the process noise here is negligible)
    const AttitudeEKF::Mat6 &after = nominal.Covariance();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            decimal expected = 0;
            for (int k = 0; k < 6; k++) {
                for (int l = 0; l < 6; l++) expected += phi[i][k] * before[k][l] * phi[j][l];
            }
            INFO("covariance entry " << i << ", " << j);
            CHECK(after[i][j] == Approx(expected).margin(2e-8));
        }
    }
}

TEST_CASE("EKF covariance after an update is expressed about the corrected attitude", "[kalman]") {
    // a large correction and very different sigmas, so that the change of reference shows
    const decimal rollSigma = DECIMAL(0.6), crossSigma = DECIMAL(0.2), baseline = DECIMAL(0.5);
    const decimal variance[3] = {rollSigma * rollSigma, crossSigma * crossSigma, crossSigma * crossSigma};
    const Quaternion identity(1, 0, 0, 0);
    AttitudeEKF ekf(crossSigma, rollSigma, DECIMAL(0.01));
    REQUIRE(ekf.Reset(identity, Vec3(0, 0, 0), baseline));
    REQUIRE(ekf.Update(RotationVectorToQuaternion(Vec3(0, 0, DECIMAL(0.6))), DECIMAL(16.27)));

    // Straight after a reset each axis is independent, with as much prior as measurement variance,
    // so the update goes halfway to the measurement...
    const Vec3 correction(0, 0, DECIMAL(0.3));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), RotationVectorToQuaternion(correction)) < DECIMAL(1e-5));
    // ...and leaves this covariance, as seen from the attitude before the correction
    AttitudeEKF::Mat6 aboutOld{};
    for (int i = 0; i < 3; i++) {
        aboutOld[i][i] = variance[i] / 2;
        aboutOld[i][i + 3] = aboutOld[i + 3][i] = variance[i] / (2 * baseline);
        aboutOld[i + 3][i + 3] = 3 * variance[i] / (2 * baseline * baseline);
    }

    // change of error coordinates from the old attitude to the corrected one, by central differences
    AttitudeEKF::Mat6 change{};
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const decimal epsilon = DECIMAL(1e-2);
    for (int j = 0; j < 3; j++) {
        const Quaternion plus = RotationVectorToQuaternion(correction + basis[j] * epsilon);
        const Quaternion minus = RotationVectorToQuaternion(correction - basis[j] * epsilon);
        const Vec3 column = QuaternionToRotationVector(plus * minus.Conjugate()) * (1 / (2 * epsilon));
        change[0][j] = column.x;
        change[1][j] = column.y;
        change[2][j] = column.z;
        change[j + 3][j + 3] = 1;
    }

    const AttitudeEKF::Mat6 &after = ekf.Covariance();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            decimal expected = 0;
            for (int k = 0; k < 6; k++) {
                for (int l = 0; l < 6; l++) expected += change[i][k] * aboutOld[k][l] * change[j][l];
            }
            INFO("covariance entry " << i << ", " << j);
            CHECK(after[i][j] == Approx(expected).margin(1e-4));
        }
    }
    // roll and cross-boresight uncertainty differ, so the turn about z couples x and y
    CHECK(after[0][1] == Approx(0.0235).margin(5e-4));
}

TEST_CASE("Rotation vector conversions match known values at every scale", "[kalman]") {
    // a quarter turn about z
    const Quaternion quarter = RotationVectorToQuaternion(Vec3(0, 0, DECIMAL(M_PI) / 2));
    CHECK(quarter.real == Approx(std::sqrt(0.5)).margin(1e-6));
    CHECK(quarter.i == Approx(0.0).margin(1e-6));
    CHECK(quarter.j == Approx(0.0).margin(1e-6));
    CHECK(quarter.k == Approx(std::sqrt(0.5)).margin(1e-6));
    // one radian about x
    const Vec3 one = QuaternionToRotationVector(Quaternion(DECIMAL_COS(DECIMAL(0.5)), DECIMAL_SIN(DECIMAL(0.5)), 0, 0));
    CHECK(one.x == Approx(1.0).margin(1e-6));
    CHECK(one.y == Approx(0.0).margin(1e-6));
    CHECK(one.z == Approx(0.0).margin(1e-6));

    const Vec3 direction(DECIMAL(0.6), DECIMAL(-0.48), DECIMAL(0.64));
    const decimal magnitudes[] = {DECIMAL(1e-13), DECIMAL(1e-9), DECIMAL(1e-5), DECIMAL(1e-3), DECIMAL(0.01),
                                  DECIMAL(0.05), DECIMAL(0.15), DECIMAL(0.5), DECIMAL(1.0), DECIMAL(3.0)};
    for (decimal magnitude : magnitudes) {
        INFO("rotation of " << magnitude << " rad");
        const Vec3 v = direction * magnitude;
        const Quaternion q = RotationVectorToQuaternion(v);
        CHECK(NormSq(q) == Approx(1.0).margin(1e-6));
        // the vector part is the axis times sin(angle / 2)
        const decimal sinHalf = DECIMAL_SIN(magnitude / 2);
        CHECK(q.i == Approx(direction.x * sinHalf).epsilon(1e-5));
        CHECK(q.j == Approx(direction.y * sinHalf).epsilon(1e-5));
        CHECK(q.k == Approx(direction.z * sinHalf).epsilon(1e-5));
        const Vec3 back = QuaternionToRotationVector(q);
        CHECK(back.x == Approx(v.x).epsilon(1e-5));
        CHECK(back.y == Approx(v.y).epsilon(1e-5));
        CHECK(back.z == Approx(v.z).epsilon(1e-5));
        const Vec3 negated = QuaternionToRotationVector(Quaternion(-q.real, -q.i, -q.j, -q.k));
        CHECK(negated.x == Approx(v.x).epsilon(1e-5));
        CHECK(negated.y == Approx(v.y).epsilon(1e-5));
        CHECK(negated.z == Approx(v.z).epsilon(1e-5));
    }

    // exactly half a turn: q and -q give opposite vectors, each a rotation of pi about the same line
    const Vec3 half = QuaternionToRotationVector(Quaternion(0, 1, 0, 0));
    const Vec3 halfNegated = QuaternionToRotationVector(Quaternion(0, -1, 0, 0));
    CHECK(half.x == Approx(M_PI).epsilon(1e-6));
    CHECK(halfNegated.x == Approx(-M_PI).epsilon(1e-6));
    CHECK(half.y == 0);
    CHECK(half.z == 0);
}

TEST_CASE("EKF Reset from two fixes with no turn between them stores a unit attitude and R, R/b, 2R/b^2 + s*b/3", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec, baseline = DECIMAL(0.25);
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.3));
    // within kUnitTolerance of unit, but not unit
    const decimal scale = DECIMAL(1.0002);
    const Quaternion scaled(start.real * scale, start.i * scale, start.j * scale, start.k * scale);
    const Vec3 rate(0, 0, 0);
    // the rate is the average over the baseline; angular acceleration makes the rate at its end differ
    const decimal accelerationSigma = DECIMAL(5e-4);
    const decimal wander = accelerationSigma * accelerationSigma * baseline / 3;
    // it is a real part of the rate variance: more than the cross-boresight fix errors' share
    CHECK(wander > 2 * crossSigma * crossSigma / (baseline * baseline));

    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
    REQUIRE(ekf.Reset(scaled, rate, baseline));
    CHECK(NormSq(ekf.AttitudeQuaternion()) == Approx(1.0).margin(1e-6));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), start) < DECIMAL(1e-6));
    CHECK((ekf.Rate() - rate).Magnitude() == 0);

    const decimal variance[3] = {rollSigma * rollSigma, crossSigma * crossSigma, crossSigma * crossSigma};
    const AttitudeEKF::Mat6 &covariance = ekf.Covariance();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            INFO("covariance entry " << i << ", " << j);
            if (i % 3 != j % 3) {
                CHECK(covariance[i][j] == 0);
            } else if (i == j && i >= 3) {
                CHECK(covariance[i][j] == Approx(2 * variance[i % 3] / (baseline * baseline) + wander).epsilon(1e-5));
            } else if (i < 3 && j < 3) {
                CHECK(covariance[i][j] == Approx(variance[i]).epsilon(1e-5));
            } else if (i >= 3 && j >= 3) {
                CHECK(covariance[i][j] == Approx(2 * variance[i % 3] / (baseline * baseline) + wander).epsilon(1e-5));
            } else {
                CHECK(covariance[i][j] == Approx(variance[i % 3] / baseline).epsilon(1e-5));
            }
        }
    }
}

TEST_CASE("EKF process noise is that of white angular acceleration", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const decimal accelerationSigma = DECIMAL(1e-3), baseline = DECIMAL(0.5), dt = DECIMAL(0.3);
    const decimal a = accelerationSigma * accelerationSigma;
    const decimal variance[3] = {rollSigma * rollSigma, crossSigma * crossSigma, crossSigma * crossSigma};
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.3));

    // At zero rate the transition is exactly [I, dt * I; 0, I], which leaves the process noise in plain view
    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
    REQUIRE(ekf.Reset(start, Vec3(0, 0, 0), baseline));
    REQUIRE(ekf.Predict(dt));
    const AttitudeEKF::Mat6 &covariance = ekf.Covariance();
    for (int i = 0; i < 3; i++) {
        const decimal attitude = variance[i], cross = variance[i] / baseline;
        const decimal rate = 2 * variance[i] / (baseline * baseline) + a * baseline / 3;
        INFO("axis " << i);
        // the process noise is a large part of each: about 9e-9, 4.5e-8 and 3e-7
        CHECK(covariance[i][i] == Approx(attitude + 2 * dt * cross + dt * dt * rate + a * dt * dt * dt / 3).epsilon(1e-4));
        CHECK(covariance[i][i + 3] == Approx(cross + dt * rate + a * dt * dt / 2).epsilon(1e-4));
        CHECK(covariance[i + 3][i] == Approx(cross + dt * rate + a * dt * dt / 2).epsilon(1e-4));
        CHECK(covariance[i + 3][i + 3] == Approx(rate + a * dt).epsilon(1e-4));
    }
    CHECK(covariance[0][1] == 0);
    CHECK(covariance[1][5] == 0);

    // two half steps add the same noise as one whole step
    AttitudeEKF halves(crossSigma, rollSigma, accelerationSigma);
    REQUIRE(halves.Reset(start, Vec3(0, 0, 0), baseline));
    REQUIRE(halves.Predict(dt / 2));
    REQUIRE(halves.Predict(dt / 2));
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            INFO("covariance entry " << i << ", " << j);
            CHECK(halves.Covariance()[i][j] == Approx(covariance[i][j]).epsilon(1e-4));
        }
    }
}

TEST_CASE("EKF covariance prediction matches finite differences for several rates and steps", "[kalman]") {
    struct Scenario { Vec3 rate; decimal dt; };
    const Scenario scenarios[] = {
        {Vec3(DECIMAL(0.3), DECIMAL(-0.2), DECIMAL(1.5)), DECIMAL(0.25)},  // dt != 1: w and w * dt differ
        {Vec3(DECIMAL(0.3), DECIMAL(-0.2), DECIMAL(0.5)), DECIMAL(0.1)},   // a 3.5 degree step
        {Vec3(DECIMAL(0.3), DECIMAL(-0.2), DECIMAL(1.5)), DECIMAL(2.0)},   // more than half a turn
    };
    const Quaternion start = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));
    const decimal baseline = DECIMAL(0.5), accelerationSigma = DECIMAL(1e-5), epsilon = DECIMAL(1e-2);
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    for (const Scenario &scenario : scenarios) {
        const Vec3 rate = scenario.rate;
        const decimal dt = scenario.dt;
        auto predicted = [&](const Quaternion &attitude, const Vec3 &startRate) {
            AttitudeEKF ekf(DECIMAL(1e-3), DECIMAL(3e-3), accelerationSigma);
            REQUIRE(ekf.Reset(attitude, startRate, baseline));
            REQUIRE(ekf.Predict(dt));
            return ekf;
        };
        const AttitudeEKF nominal = predicted(start, rate);
        AttitudeEKF atReset(DECIMAL(1e-3), DECIMAL(3e-3), accelerationSigma);
        REQUIRE(atReset.Reset(start, rate, baseline));
        const AttitudeEKF::Mat6 before = atReset.Covariance();

        AttitudeEKF::Mat6 phi{};
        for (int j = 0; j < 6; j++) {
            const Quaternion plusAttitude = j < 3 ? RotationVectorToQuaternion(basis[j] * epsilon) * start : start;
            const Quaternion minusAttitude = j < 3 ? RotationVectorToQuaternion(basis[j] * -epsilon) * start : start;
            const Vec3 plusRate = j < 3 ? rate : rate + basis[j - 3] * epsilon;
            const Vec3 minusRate = j < 3 ? rate : rate - basis[j - 3] * epsilon;
            const AttitudeEKF plus = predicted(plusAttitude, plusRate);
            const AttitudeEKF minus = predicted(minusAttitude, minusRate);
            const Vec3 attitudeChange = QuaternionToRotationVector(
                plus.AttitudeQuaternion() * minus.AttitudeQuaternion().Conjugate()) * (1 / (2 * epsilon));
            const Vec3 rateChange = (plus.Rate() - minus.Rate()) * (1 / (2 * epsilon));
            const decimal column[6] = {attitudeChange.x, attitudeChange.y, attitudeChange.z,
                                       rateChange.x, rateChange.y, rateChange.z};
            for (int i = 0; i < 6; i++) phi[i][j] = column[i];
        }

        const AttitudeEKF::Mat6 &after = nominal.Covariance();
        decimal largest = 0;
        for (int i = 0; i < 6; i++) largest = std::max(largest, after[i][i]);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                decimal expected = 0;
                for (int k = 0; k < 6; k++) {
                    for (int l = 0; l < 6; l++) expected += phi[i][k] * before[k][l] * phi[j][l];
                }
                INFO("dt " << dt << ", covariance entry " << i << ", " << j);
                CHECK(after[i][j] == Approx(expected).margin(largest * DECIMAL(3e-4)));
            }
        }
    }
}

TEST_CASE("EKF update with a correlated attitude covariance matches the textbook update", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const double variance[3] = {static_cast<double>(rollSigma * rollSigma), static_cast<double>(crossSigma * crossSigma),
                                static_cast<double>(crossSigma * crossSigma)};
    const Quaternion start = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));

    AttitudeEKF ekf(crossSigma, rollSigma, DECIMAL(1e-5));
    // a long baseline keeps the prior about as large as the measurement noise, so the gain is far from 0 and 1;
    // the turn of 0.8 rad about a skew axis mixes roll into the other axes
    REQUIRE(ekf.Reset(start, Vec3(DECIMAL(0.6), DECIMAL(-1.0), DECIMAL(1.2)), DECIMAL(0.9)));
    REQUIRE(ekf.Predict(DECIMAL(0.5)));
    const AttitudeEKF::Mat6 prior = ekf.Covariance();
    const Quaternion priorAttitude = ekf.AttitudeQuaternion();
    const Vec3 priorRate = ekf.Rate();
    // the attitude block really is correlated
    CHECK(std::abs(prior[0][1]) > DECIMAL(0.1) * prior[1][1]);
    CHECK(std::abs(prior[1][2]) > DECIMAL(0.1) * prior[1][1]);

    const Vec3 innovation(20 * kArcsec, -8 * kArcsec, 5 * kArcsec);
    REQUIRE(ekf.Update(RotationVectorToQuaternion(innovation) * priorAttitude, DECIMAL(16.27)));

    // reference in double precision: S = P_theta + R, inverted by cofactors
    double s[3][3], sInverse[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) s[i][j] = static_cast<double>(prior[i][j]) + (i == j ? variance[i] : 0.0);
    }
    const double determinant = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1])
        - s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0])
        + s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            const int r0 = (j + 1) % 3, r1 = (j + 2) % 3, c0 = (i + 1) % 3, c1 = (i + 2) % 3;
            sInverse[i][j] = (s[r0][c0] * s[r1][c1] - s[r0][c1] * s[r1][c0]) / determinant;
        }
    }
    // K = P H^T S^-1; correction = K * innovation; P+ = P - K H P
    const double innovationArray[3] = {static_cast<double>(innovation.x), static_cast<double>(innovation.y),
                                       static_cast<double>(innovation.z)};
    double gain[6][3] = {}, correction[6] = {};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) gain[i][j] += static_cast<double>(prior[i][k]) * sInverse[k][j];
            correction[i] += gain[i][j] * innovationArray[j];
        }
    }
    const Quaternion expectedAttitude = RotationVectorToQuaternion(
        Vec3(static_cast<decimal>(correction[0]), static_cast<decimal>(correction[1]), static_cast<decimal>(correction[2])))
        * priorAttitude;
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), expectedAttitude) < kArcsec / 10);
    const Vec3 rateChange = ekf.Rate() - priorRate;
    // the rate is near 1 rad/s here, so in float mode its change is only resolved to about 1e-7 rad/s
    const double rateResolution = 8 * static_cast<double>(std::numeric_limits<decimal>::epsilon() * priorRate.Magnitude());
    CHECK(static_cast<double>(rateChange.x) == Approx(correction[3]).epsilon(1e-3).margin(rateResolution));
    CHECK(static_cast<double>(rateChange.y) == Approx(correction[4]).epsilon(1e-3).margin(rateResolution));
    CHECK(static_cast<double>(rateChange.z) == Approx(correction[5]).epsilon(1e-3).margin(rateResolution));
    // the corrections are not trivially small: several arcseconds, and arcseconds per second
    CHECK(std::abs(correction[0]) > static_cast<double>(kArcsec));
    CHECK(std::abs(correction[4]) > static_cast<double>(kArcsec));

    const AttitudeEKF::Mat6 &posterior = ekf.Covariance();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            double expected = static_cast<double>(prior[i][j]);
            for (int k = 0; k < 3; k++) expected -= gain[i][k] * static_cast<double>(prior[k][j]);
            INFO("covariance entry " << i << ", " << j);
            // the re-centring Jacobian differs from the identity by 1e-4 here, well inside the tolerance
            CHECK(static_cast<double>(posterior[i][j]) == Approx(expected).epsilon(2e-3).margin(2e-3 * variance[1]));
        }
    }
}

TEST_CASE("EKF gate sits at the stated Mahalanobis distance on each axis", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec, gate = DECIMAL(16.27);
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.3));
    // straight after a Reset the innovation variance on each axis is twice the measurement variance
    const decimal rollLimit = DECIMAL_SQRT(gate * 2) * rollSigma, crossLimit = DECIMAL_SQRT(gate * 2) * crossSigma;
    auto accepts = [&](const Vec3 &offset) {
        AttitudeEKF ekf(crossSigma, rollSigma, DECIMAL(0.01));
        REQUIRE(ekf.Reset(start, Vec3(0, 0, 0), DECIMAL(0.1)));
        const AttitudeEKF::Mat6 covarianceBefore = ekf.Covariance();
        const Quaternion attitudeBefore = ekf.AttitudeQuaternion();
        const bool accepted = ekf.Update(RotationVectorToQuaternion(offset) * start, gate);
        if (!accepted) {
            // a rejected measurement changes nothing at all
            CHECK(AngleBetween(ekf.AttitudeQuaternion(), attitudeBefore) == 0);
            CHECK(ekf.Rate().Magnitude() == 0);
            CHECK(ekf.Covariance() == covarianceBefore);
        }
        return accepted;
    };
    CHECK(accepts(Vec3(rollLimit * DECIMAL(0.98), 0, 0)));
    CHECK_FALSE(accepts(Vec3(rollLimit * DECIMAL(1.02), 0, 0)));
    CHECK(accepts(Vec3(0, crossLimit * DECIMAL(0.98), 0)));
    CHECK_FALSE(accepts(Vec3(0, crossLimit * DECIMAL(1.02), 0)));
    CHECK(accepts(Vec3(0, 0, -crossLimit * DECIMAL(0.98))));
    CHECK_FALSE(accepts(Vec3(0, 0, -crossLimit * DECIMAL(1.02))));
    // all three axes at once, each at 1/sqrt(3) of its limit
    const decimal third = DECIMAL_SQRT(DECIMAL(1.0) / 3);
    CHECK(accepts(Vec3(rollLimit * third * DECIMAL(0.98), crossLimit * third * DECIMAL(0.98), crossLimit * third * DECIMAL(0.98))));
    CHECK_FALSE(accepts(Vec3(rollLimit * third * DECIMAL(1.02), crossLimit * third * DECIMAL(1.02), crossLimit * third * DECIMAL(1.02))));
}

TEST_CASE("EKF state stays unit, symmetric and positive definite over a long run", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec, dt = DECIMAL(0.1);
    const Vec3 rate(DECIMAL(0.02), DECIMAL(-0.03), DECIMAL(0.01));
    Quaternion truth = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.3));
    const Quaternion step = RotationVectorToQuaternion(rate * dt);
    Gaussian gaussian(7);

    AttitudeEKF ekf(crossSigma, rollSigma, DECIMAL(1e-4));
    REQUIRE(ekf.Reset(truth, rate, dt));
    const decimal unitTolerance = 8 * std::numeric_limits<decimal>::epsilon();
    int accepted = 0;
    for (int frame = 0; frame < 20000; frame++) {
        truth = step * truth;
        REQUIRE(ekf.Predict(dt));
        const Vec3 noise(rollSigma * gaussian(), crossSigma * gaussian(), crossSigma * gaussian());
        if (ekf.Update(RotationVectorToQuaternion(noise) * truth, DECIMAL(16.27))) accepted++;

        const AttitudeEKF::Mat6 &covariance = ekf.Covariance();
        bool symmetric = true;
        const bool positive = IsPositiveDefinite(covariance, 1e-6);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < i; j++) symmetric = symmetric && covariance[i][j] == covariance[j][i];
        }
        if (!symmetric || !positive || std::abs(NormSq(ekf.AttitudeQuaternion()) - 1) > unitTolerance) {
            INFO("frame " << frame);
            REQUIRE(symmetric);
            REQUIRE(positive);
            REQUIRE(std::abs(NormSq(ekf.AttitudeQuaternion()) - 1) <= unitTolerance);
        }
    }
    // the gate rejects 0.1% of honest measurements
    CHECK(accepted > 19900);
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), truth) < 30 * kArcsec);
    // prediction alone, as when no fix arrives for a long time, also keeps the attitude unit
    // and the covariance positive definite
    for (int frame = 0; frame < 20000; frame++) {
        REQUIRE(ekf.Predict(dt));
        if (!IsPositiveDefinite(ekf.Covariance(), 1e-6)) {
            INFO("frame " << frame << " of prediction alone");
            REQUIRE(IsPositiveDefinite(ekf.Covariance(), 1e-6));
        }
    }
    CHECK(std::abs(NormSq(ekf.AttitudeQuaternion()) - 1) <= unitTolerance);
}

TEST_CASE("EKF covariance is consistent with its actual error over many noisy runs", "[kalman]") {
    // A star tracker at 10 Hz turning at 1.5 degrees per second, with 5 and 30 arcsecond fixes and a
    // randomly wandering rate. Normalized estimation error squared (NEES) of a consistent 6-state
    // filter averages 6.
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const decimal accelerationSigma = DECIMAL(2e-4), dt = DECIMAL(0.1);
    const int runs = 300, frames = 100, firstScored = 20;
    Gaussian gaussian(20261120);

    double neesSum = 0, attitudeSquared[3] = {}, rateSquared = 0;
    int samples = 0, rejected = 0;
    for (int run = 0; run < runs; run++) {
        Quaternion truth = SphericalToQuaternion(DECIMAL(0.02) * run, DECIMAL(0.3), DECIMAL(0.7));
        Vec3 rate(DECIMAL(0.015), DECIMAL(-0.02), DECIMAL(0.008));
        auto advance = [&]() {
            // white angular acceleration over dt: rate change and the attitude change it causes
            Vec3 rateStep, angleStep;
            decimal *rateParts[3] = {&rateStep.x, &rateStep.y, &rateStep.z};
            decimal *angleParts[3] = {&angleStep.x, &angleStep.y, &angleStep.z};
            for (int i = 0; i < 3; i++) {
                const decimal n1 = gaussian(), n2 = gaussian();
                *rateParts[i] = accelerationSigma * DECIMAL_SQRT(dt) * n1;
                *angleParts[i] = accelerationSigma * dt * DECIMAL_SQRT(dt) * (n1 / 2 + n2 / DECIMAL_SQRT(DECIMAL(12.0)));
            }
            truth = RotationVectorToQuaternion(rate * dt + angleStep) * truth;
            rate = rate + rateStep;
        };
        auto measure = [&]() {
            const Vec3 noise(rollSigma * gaussian(), crossSigma * gaussian(), crossSigma * gaussian());
            return RotationVectorToQuaternion(noise) * truth;
        };

        // start as documented: two fixes one frame apart
        const Quaternion first = measure();
        advance();
        const Quaternion second = measure();
        AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
        REQUIRE(ekf.Reset(second, QuaternionToRotationVector(second * first.Conjugate()) * (1 / dt), dt));

        for (int frame = 0; frame < frames; frame++) {
            advance();
            REQUIRE(ekf.Predict(dt));
            if (!ekf.Update(measure(), DECIMAL(16.27))) rejected++;
            if (frame < firstScored) continue;
            const Vec3 attitudeError = QuaternionToRotationVector(truth * ekf.AttitudeQuaternion().Conjugate());
            const Vec3 rateError = rate - ekf.Rate();
            const double error[6] = {static_cast<double>(attitudeError.x), static_cast<double>(attitudeError.y),
                                     static_cast<double>(attitudeError.z), static_cast<double>(rateError.x),
                                     static_cast<double>(rateError.y), static_cast<double>(rateError.z)};
            neesSum += Nees(ekf.Covariance(), error);
            for (int i = 0; i < 3; i++) attitudeSquared[i] += error[i] * error[i];
            rateSquared += error[3] * error[3] + error[4] * error[4] + error[5] * error[5];
            samples++;
        }
    }
    const double nees = neesSum / samples;
    const double arcsec = static_cast<double>(kArcsec);
    const double rollRms = std::sqrt(attitudeSquared[0] / samples) / arcsec;
    const double crossRms = std::sqrt((attitudeSquared[1] + attitudeSquared[2]) / (2 * samples)) / arcsec;
    const double rateRms = std::sqrt(rateSquared / (3 * samples)) / arcsec;
    CAPTURE(nees, rollRms, crossRms, rateRms, rejected);
    CHECK(nees > 5.7);
    CHECK(nees < 6.3);
    // the filter improves on the raw fixes
    CHECK(rollRms < 30 * 0.8);
    CHECK(crossRms < 5 * 0.95);
    // about 0.1% of honest measurements fall outside the gate: 30 expected of 30000
    CHECK(rejected < 80);
}

TEST_CASE("EKF Reset covariance is what two noisy fixes a turn apart imply", "[kalman]") {
    // Sigmas much larger than a star tracker's, and a 40 degree turn, so that every term matters
    const decimal crossSigma = DECIMAL(1e-3), rollSigma = DECIMAL(6e-3), baseline = DECIMAL(0.5);
    const double variance[3] = {static_cast<double>(rollSigma * rollSigma), static_cast<double>(crossSigma * crossSigma),
                                static_cast<double>(crossSigma * crossSigma)};
    const Vec3 turn(DECIMAL(0.3), DECIMAL(-0.5), DECIMAL(0.4));
    const Quaternion older = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));
    const Quaternion newer = RotationVectorToQuaternion(turn) * older;
    const Vec3 basis[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // the rate a tracker would compute from two fixes with the given errors
    auto rateFrom = [&](const Vec3 &newerError, const Vec3 &olderError) {
        const Quaternion measuredNewer = RotationVectorToQuaternion(newerError) * newer;
        const Quaternion measuredOlder = RotationVectorToQuaternion(olderError) * older;
        return QuaternionToRotationVector(measuredNewer * measuredOlder.Conjugate()) * (1 / baseline);
    };
    const Vec3 rate = rateFrom(Vec3(0, 0, 0), Vec3(0, 0, 0));

    const decimal accelerationSigma = DECIMAL(0.02);
    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
    REQUIRE(ekf.Reset(newer, rate, baseline));
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), newer) < DECIMAL(1e-6));
    CHECK((ekf.Rate() - rate).Magnitude() == 0);

    // sensitivity of the rate to each fix error, by central differences; the attitude error is the
    // newer fix's error itself
    const decimal epsilon = DECIMAL(1e-2);
    double sensitivity[6][6] = {};  // rows: attitude, rate; columns: newer fix error, older fix error
    for (int j = 0; j < 6; j++) {
        const Vec3 step = basis[j % 3] * epsilon;
        const Vec3 change = j < 3
            ? (rateFrom(step, Vec3(0, 0, 0)) - rateFrom(step * -1, Vec3(0, 0, 0))) * (1 / (2 * epsilon))
            : (rateFrom(Vec3(0, 0, 0), step) - rateFrom(Vec3(0, 0, 0), step * -1)) * (1 / (2 * epsilon));
        sensitivity[3][j] = static_cast<double>(change.x);
        sensitivity[4][j] = static_cast<double>(change.y);
        sensitivity[5][j] = static_cast<double>(change.z);
        if (j < 3) sensitivity[j][j] = 1;
    }
    double expected[6][6] = {};
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) expected[i][j] += sensitivity[i][k] * variance[k % 3] * sensitivity[j][k];
        }
    }
    // The rate is the average over the baseline, used as the rate at the newer fix. An impulse of
    // angular acceleration a time t after the older fix, with the rate at the newer fix held at
    // `rate`, means the rate before it was (rate - impulse); the rate two exact fixes would give
    // is then off by this much:
    auto rateErrorFrom = [&](const Vec3 &impulse, decimal t) {
        const Quaternion turned = RotationVectorToQuaternion(rate * (baseline - t))
            * RotationVectorToQuaternion((rate - impulse) * t);
        return QuaternionToRotationVector(turned) * (1 / baseline) - rate;
    };
    // white acceleration is a train of independent impulses, each of variance density * slice
    const int slices = 400;
    const double slice = static_cast<double>(baseline) / slices;
    const double density = static_cast<double>(accelerationSigma * accelerationSigma);
    double wander[3][3] = {};
    for (int n = 0; n < slices; n++) {
        const decimal t = static_cast<decimal>((n + 0.5) * slice);
        double impulseSensitivity[3][3];
        for (int j = 0; j < 3; j++) {
            const Vec3 change = (rateErrorFrom(basis[j] * epsilon, t) - rateErrorFrom(basis[j] * -epsilon, t))
                * (1 / (2 * epsilon));
            impulseSensitivity[0][j] = static_cast<double>(change.x);
            impulseSensitivity[1][j] = static_cast<double>(change.y);
            impulseSensitivity[2][j] = static_cast<double>(change.z);
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) wander[i][j] += density * slice * impulseSensitivity[i][k] * impulseSensitivity[j][k];
            }
        }
    }
    // both parts matter here: the wander is between a tenth and ten times the fix errors' share
    CHECK(wander[1][1] > 0.1 * expected[4][4]);
    CHECK(wander[1][1] < 10 * expected[4][4]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) expected[i + 3][j + 3] += wander[i][j];
    }

    const AttitudeEKF::Mat6 &covariance = ekf.Covariance();
    CHECK(IsPositiveDefinite(covariance, 1e-6));
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            INFO("covariance entry " << i << ", " << j);
            CHECK(static_cast<double>(covariance[i][j])
                  == Approx(expected[i][j]).margin(1e-3 * std::sqrt(expected[i][i] * expected[j][j])));
        }
    }
    // the turn matters: the rate variance about z is far from the 2R/b^2 of two fixes with no turn
    CHECK(expected[5][5] - wander[2][2] > 3 * 2 * variance[2] / static_cast<double>(baseline * baseline));
}

TEST_CASE("EKF prediction covariance matches numerical integration through a turn", "[kalman]") {
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const Quaternion start = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));
    decimal baseline = DECIMAL(0.5);
    Vec3 rate;
    decimal dt, accelerationSigma;
    bool noiseDominates;
    double tolerance = 1e-3;  // of each entry's scale
    SECTION("a 115 degree step, as after a long gap") {
        // the process noise is a three-point quadrature, so it is approximate for a turn this large
        rate = Vec3(DECIMAL(0.012), DECIMAL(-0.010), DECIMAL(0.013));
        dt = DECIMAL(100.0);
        accelerationSigma = DECIMAL(2e-4);
        noiseDominates = true;
        tolerance = 1e-3;
    }
    SECTION("a 3.5 degree step with strong process noise") {
        // a rate far beyond a star tracker's, so that the turn's effect on the process noise is in plain view
        rate = Vec3(DECIMAL(0.3), DECIMAL(-0.4), DECIMAL(0.35));
        dt = DECIMAL(0.1);
        accelerationSigma = DECIMAL(0.05);
        noiseDominates = true;
        baseline = DECIMAL(0.02);  // short, so that the rate uncertainty Reset leaves is small beside the noise
    }
    SECTION("a 0.45 degree step with almost no process noise") {
        // a step small enough for the Jacobian's series, over a time long enough for it to matter
        rate = Vec3(DECIMAL(0.004), DECIMAL(-0.005), DECIMAL(0.0045));
        dt = DECIMAL(1.0);
        accelerationSigma = DECIMAL(1e-7);
        noiseDominates = false;
    }

    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
    REQUIRE(ekf.Reset(start, rate, baseline));
    const AttitudeEKF::Mat6 before = ekf.Covariance();
    REQUIRE(ekf.Predict(dt));
    const AttitudeEKF::Mat6 &after = ekf.Covariance();

    // The error dynamics are d(theta)/dt = [w]x * theta + dw, d(dw)/dt = white noise. With
    // E(t) the rotation of the turn over t and B(t) the integral of E from 0 to t:
    //   transition = [E(dt), B(dt); 0, I],
    //   process noise = density * integral over t of [B B^T, B; B^T, I].
    // Both integrals by the trapezoid rule.
    const int steps = 2000;
    const double h = static_cast<double>(dt) / steps;
    const double density = static_cast<double>(accelerationSigma * accelerationSigma);
    double e[3][3], b[3][3] = {}, noise[6][6] = {};
    RotationMatrix(Quaternion(1, 0, 0, 0), e);
    auto accumulate = [&](double weight) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) noise[i][j] += weight * b[i][k] * b[j][k];
                noise[i][j + 3] += weight * b[i][j];
                noise[j + 3][i] += weight * b[i][j];
            }
            noise[i + 3][i + 3] += weight;
        }
    };
    accumulate(h / 2);
    for (int step = 1; step <= steps; step++) {
        double next[3][3];
        RotationMatrix(RotationVectorToQuaternion(rate * static_cast<decimal>(step * h)), next);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                b[i][j] += h / 2 * (e[i][j] + next[i][j]);
                e[i][j] = next[i][j];
            }
        }
        accumulate(step == steps ? h / 2 : h);
    }
    double transition[6][6] = {};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            transition[i][j] = e[i][j];
            transition[i][j + 3] = b[i][j];
        }
        transition[i + 3][i + 3] = 1;
    }
    double expected[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            expected[i][j] = density * noise[i][j];
            for (int k = 0; k < 6; k++) {
                for (int l = 0; l < 6; l++) {
                    expected[i][j] += transition[i][k] * static_cast<double>(before[k][l]) * transition[j][l];
                }
            }
        }
    }
    CHECK(IsPositiveDefinite(after, 1e-6));
    double worst = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            INFO("covariance entry " << i << ", " << j);
            const double scale = std::sqrt(expected[i][i] * expected[j][j]);
            worst = std::max(worst, std::abs(static_cast<double>(after[i][j]) - expected[i][j]) / scale);
            CHECK(static_cast<double>(after[i][j]) == Approx(expected[i][j]).margin(tolerance * scale));
        }
    }
    CAPTURE(worst);
    // each scenario shows what it is there for: the process noise dominates, or it is negligible
    const double noiseShare = density * noise[4][4] / expected[4][4];
    CHECK((noiseDominates ? noiseShare > 0.9 : noiseShare < 0.01));
    // and the turn shows: the attitude-rate block is not symmetric
    CHECK(std::abs(expected[0][4] - expected[1][3]) > 5e-3 * std::sqrt(expected[0][0] * expected[4][4]));
}

TEST_CASE("EKF covariance stays positive definite through one long step and the update after it", "[kalman]") {
    // The settings of the consistency test, turning at 1.15 degrees per second, with no fix for
    // a while: one Predict call covers the whole gap.
    const decimal gaps[] = {DECIMAL(1.0), DECIMAL(10.0), DECIMAL(50.0), DECIMAL(87.0), DECIMAL(100.0),
                            DECIMAL(157.0), DECIMAL(400.0), DECIMAL(3000.0)};
    for (decimal gap : gaps) {
        INFO("gap of " << gap << " s, a turn of " << gap * DECIMAL(0.02) << " rad");
        AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(2e-4));
        REQUIRE(ekf.Reset(Quaternion(1, 0, 0, 0), Vec3(0, 0, DECIMAL(0.02)), DECIMAL(1.0)));
        REQUIRE(ekf.Predict(gap));
        CHECK(IsPositiveDefinite(ekf.Covariance(), 1e-6));
        // the measurement agrees exactly with the prediction
        REQUIRE(ekf.Update(ekf.AttitudeQuaternion(), DECIMAL(16.27)));
        CHECK(IsPositiveDefinite(ekf.Covariance(), 1e-6));
        for (int i = 0; i < 6; i++) CHECK(ekf.Covariance()[i][i] > 0);
    }
}

TEST_CASE("EKF prediction over a gap is nearly the same in one step as in many", "[kalman]") {
    struct Gap { decimal seconds; double tolerance; };
    // turning at 1.16 degrees per second: turns of 0.2, 0.5, 1, 2, 6.4 and 20 radians
    const Gap gaps[] = {{DECIMAL(10.0), 1e-3}, {DECIMAL(25.0), 1e-3}, {DECIMAL(50.0), 1e-3}, {DECIMAL(100.0), 1e-3},
                        {DECIMAL(314.0), 1e-3}, {DECIMAL(1000.0), 1e-3}};
    const Vec3 rate(DECIMAL(0.012), DECIMAL(-0.010), DECIMAL(0.013));
    const Quaternion start = SphericalToQuaternion(DECIMAL(0.4), DECIMAL(-0.3), DECIMAL(1.1));
    for (const Gap &gap : gaps) {
        INFO("gap of " << gap.seconds << " s");
        AttitudeEKF whole(5 * kArcsec, 30 * kArcsec, DECIMAL(2e-4));
        AttitudeEKF stepped(5 * kArcsec, 30 * kArcsec, DECIMAL(2e-4));
        REQUIRE(whole.Reset(start, rate, DECIMAL(1.0)));
        REQUIRE(stepped.Reset(start, rate, DECIMAL(1.0)));
        REQUIRE(whole.Predict(gap.seconds));
        const int steps = 1000;
        for (int step = 0; step < steps; step++) REQUIRE(stepped.Predict(gap.seconds / steps));

        CHECK(AngleBetween(whole.AttitudeQuaternion(), stepped.AttitudeQuaternion()) < DECIMAL(1e-3));
        CHECK(IsPositiveDefinite(whole.Covariance(), 1e-6));
        CHECK(IsPositiveDefinite(stepped.Covariance(), 1e-6));
        double worst = 0;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                INFO("covariance entry " << i << ", " << j);
                const double scale = std::sqrt(static_cast<double>(stepped.Covariance()[i][i])
                                               * static_cast<double>(stepped.Covariance()[j][j]));
                const double difference = static_cast<double>(whole.Covariance()[i][j])
                    - static_cast<double>(stepped.Covariance()[i][j]);
                worst = std::max(worst, std::abs(difference) / scale);
                CHECK(std::abs(difference) <= gap.tolerance * scale);
            }
        }
        CAPTURE(worst);
    }
}

TEST_CASE("EKF covariance straight after Reset is consistent with its actual error", "[kalman]") {
    // Two fixes a second apart on a spacecraft whose rate wanders in between. The rate they give
    // is the average over that second, not the rate at its end, and the covariance has to say so:
    // here that is eleven times the uncertainty the fix errors alone put on the cross-boresight rate.
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const decimal accelerationSigma = DECIMAL(2e-4), dt = DECIMAL(0.1);
    const int runs = 2000, framesBetweenFixes = 10;
    Gaussian gaussian(20261121);

    double neesSum = 0;
    int rejectedAtFirstGate = 0;
    for (int run = 0; run < runs; run++) {
        NoisyTruth truth(SphericalToQuaternion(DECIMAL(0.003) * run, DECIMAL(0.3), DECIMAL(0.7)),
                         Vec3(DECIMAL(0.015), DECIMAL(-0.02), DECIMAL(0.008)),
                         crossSigma, rollSigma, accelerationSigma, &gaussian);
        const Quaternion first = truth.Measure();
        for (int frame = 0; frame < framesBetweenFixes; frame++) truth.Advance(dt);
        const Quaternion second = truth.Measure();
        const decimal baseline = dt * framesBetweenFixes;
        AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
        REQUIRE(ekf.Reset(second, QuaternionToRotationVector(second * first.Conjugate()) * (1 / baseline), baseline));

        const Vec3 attitudeError = QuaternionToRotationVector(truth.attitude * ekf.AttitudeQuaternion().Conjugate());
        const Vec3 rateError = truth.rate - ekf.Rate();
        const double error[6] = {static_cast<double>(attitudeError.x), static_cast<double>(attitudeError.y),
                                 static_cast<double>(attitudeError.z), static_cast<double>(rateError.x),
                                 static_cast<double>(rateError.y), static_cast<double>(rateError.z)};
        neesSum += Nees(ekf.Covariance(), error);

        // and the next fix, a second later, passes the gate as often as it should
        for (int frame = 0; frame < framesBetweenFixes; frame++) truth.Advance(dt);
        REQUIRE(ekf.Predict(baseline));
        if (!ekf.Update(truth.Measure(), DECIMAL(16.27))) rejectedAtFirstGate++;
    }
    const double nees = neesSum / runs;
    CAPTURE(nees, rejectedAtFirstGate);
    // the mean of 2000 samples of a chi-square with 6 degrees of freedom is 6 give or take 0.08
    CHECK(nees > 5.6);
    CHECK(nees < 6.4);
    // 0.1% of honest measurements fall outside the gate: 2 expected of 2000
    CHECK(rejectedAtFirstGate <= 10);
}

TEST_CASE("EKF prediction through whole turns matches the exact process noise, in one step or many", "[kalman]") {
    // Turning about z at 0.02 rad/s for exactly two turns. Over whole turns a transverse rate
    // error integrates to no transverse attitude error at all, and the process noise on a
    // transverse attitude axis is exactly 2 * density * time / rate^2.
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;
    const decimal accelerationSigma = DECIMAL(2e-6), rateZ = DECIMAL(0.02);
    const decimal gap = 4 * DECIMAL_M_PI / rateZ;  // 628 s
    const double noise = 2 * static_cast<double>(accelerationSigma * accelerationSigma) * static_cast<double>(gap)
        / static_cast<double>(rateZ * rateZ);
    const double expectedRoll = static_cast<double>(rollSigma * rollSigma) + noise;
    const double expectedCross = static_cast<double>(crossSigma * crossSigma) + noise;
    // the noise is nearly all of it, so this is a test of the noise
    CHECK(noise > 0.99 * expectedRoll);

    const Quaternion start(1, 0, 0, 0);
    const int steps = 1000;
    auto predicted = [&](bool inOneStep) {
        AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
        REQUIRE(ekf.Reset(start, Vec3(0, 0, rateZ), DECIMAL(1.0)));
        if (inOneStep) {
            REQUIRE(ekf.Predict(gap));
        } else {
            for (int step = 0; step < steps; step++) REQUIRE(ekf.Predict(gap / steps));
        }
        return ekf;
    };
    for (bool inOneStep : {true, false}) {
        INFO((inOneStep ? "one step" : "a thousand steps"));
        const AttitudeEKF ekf = predicted(inOneStep);
        CHECK(static_cast<double>(ekf.Covariance()[0][0]) == Approx(expectedRoll).epsilon(2e-3));
        CHECK(static_cast<double>(ekf.Covariance()[1][1]) == Approx(expectedCross).epsilon(2e-3));

        // So the same fix gets the same verdict either way. About x the innovation variance is
        // expectedRoll plus the roll measurement variance: 0.0123 rad is 12.0 against the gate
        // of 16.27, and 0.0150 rad is 17.8.
        AttitudeEKF inside = ekf, outside = ekf;
        CHECK(inside.Update(RotationVectorToQuaternion(Vec3(DECIMAL(0.0123), 0, 0)) * ekf.AttitudeQuaternion(),
                            DECIMAL(16.27)));
        CHECK_FALSE(outside.Update(RotationVectorToQuaternion(Vec3(DECIMAL(0.0150), 0, 0)) * ekf.AttitudeQuaternion(),
                                   DECIMAL(16.27)));
    }
}

TEST_CASE("EKF Predict takes a step of up to 64 radians of turn and refuses a longer one", "[kalman]") {
    const Quaternion start = SphericalToQuaternion(DECIMAL(1.0), DECIMAL(0.5), DECIMAL(0.3));
    const Vec3 rate(0, 0, DECIMAL(0.02));
    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(2e-4));
    REQUIRE(ekf.Reset(start, rate, DECIMAL(1.0)));
    const AttitudeEKF::Mat6 before = ekf.Covariance();

    CHECK_FALSE(ekf.Predict(DECIMAL(3300.0)));  // 66 radians
    CHECK(AngleBetween(ekf.AttitudeQuaternion(), start) < DECIMAL(1e-6));
    CHECK((ekf.Rate() - rate).Magnitude() == 0);
    CHECK(ekf.Covariance() == before);

    // the same gap in two steps
    CHECK(ekf.Predict(DECIMAL(3100.0)));  // 62 radians
    CHECK(ekf.Predict(DECIMAL(200.0)));
    CHECK(IsPositiveDefinite(ekf.Covariance(), 1e-6));
}

TEST_CASE("EKF never commits a covariance that is not positive definite", "[kalman]") {
    const Quaternion identity(1, 0, 0, 0);
    const decimal crossSigma = 5 * kArcsec, rollSigma = 30 * kArcsec;

    // Nothing moves, and no fix arrives for 273 s after two fixes 0.1 s apart, with little or no
    // process noise. The attitude and rate errors are then within a few parts in a hundred
    // million of perfectly correlated, which a float covariance cannot hold. So in float mode
    // Predict refuses; in double it goes through and the covariance is sound.
    for (decimal accelerationSigma : {DECIMAL(1e-8), DECIMAL(0.0)}) {
        INFO("acceleration sigma " << accelerationSigma);
        AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
        REQUIRE(ekf.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
        const AttitudeEKF::Mat6 before = ekf.Covariance();
        const bool predicted = ekf.Predict(DECIMAL(273.0));
#ifdef LOST_FLOAT_MODE
        CHECK_FALSE(predicted);
#else
        CHECK(predicted);
#endif
        if (predicted) {
            CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
            REQUIRE(ekf.Update(identity, DECIMAL(16.27)));
            CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
            for (int i = 0; i < 6; i++) CHECK(ekf.Covariance()[i][i] > 0);
        } else {
            // refused: nothing changed
            CHECK(ekf.Covariance() == before);
            CHECK(AngleBetween(ekf.AttitudeQuaternion(), identity) == 0);
            // a shorter step still works, and so does starting again
            CHECK(ekf.Predict(DECIMAL(1.0)));
            CHECK(ekf.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
            CHECK(ekf.Update(identity, DECIMAL(16.27)));
            CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
        }
    }

    // Whatever the baseline, the gap, the process noise and the rate: either the call is refused
    // and changes nothing, or what it leaves is positive definite.
    int committed = 0, refused = 0;
    for (decimal baseline : {DECIMAL(0.01), DECIMAL(0.1), DECIMAL(1.0)}) {
        for (decimal gap : {DECIMAL(1.0), DECIMAL(10.0), DECIMAL(100.0), DECIMAL(1e3), DECIMAL(1e4), DECIMAL(1e5),
                            DECIMAL(1e6), DECIMAL(1e8)}) {
            for (decimal accelerationSigma : {DECIMAL(0.0), DECIMAL(1e-8), DECIMAL(1e-5), DECIMAL(1e-3)}) {
                for (decimal rateZ : {DECIMAL(0.0), DECIMAL(1e-5), DECIMAL(0.02)}) {
                    INFO("baseline " << baseline << ", gap " << gap << ", acceleration sigma " << accelerationSigma
                         << ", rate " << rateZ);
                    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
                    REQUIRE(ekf.Reset(identity, Vec3(0, 0, rateZ), baseline));
                    const AttitudeEKF::Mat6 before = ekf.Covariance();
                    if (!ekf.Predict(gap)) {
                        refused++;
                        CHECK(ekf.Covariance() == before);
                        continue;
                    }
                    committed++;
                    CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
                    const AttitudeEKF::Mat6 predictedCovariance = ekf.Covariance();
                    if (ekf.Update(ekf.AttitudeQuaternion(), DECIMAL(16.27))) {
                        CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
                    } else {
                        CHECK(ekf.Covariance() == predictedCovariance);
                    }
                }
            }
        }
    }
    CAPTURE(committed, refused);
    // both outcomes are exercised
    CHECK(committed > 100);
    CHECK(refused > 20);
}


TEST_CASE("EKF Reset refuses two fixes so far apart that a fix error could reverse the rate", "[kalman]") {
    // The spacecraft turns pi - 1e-5 radians about z between two fixes 100 s apart, and the newer
    // fix is off by 2e-5 radians, 4 arcseconds, the same way. The shorter rotation between the
    // fixes is then the other way round the circle, so the rate they give has the wrong sign.
    const decimal baseline = DECIMAL(100.0);
    const Quaternion older(1, 0, 0, 0);
    const Quaternion newer = RotationVectorToQuaternion(Vec3(0, 0, DECIMAL_M_PI + DECIMAL(1e-5))) * older;
    const Vec3 rate = QuaternionToRotationVector(newer * older.Conjugate()) * (1 / baseline);
    CHECK(rate.z < 0);  // the true rate is +0.0314 rad/s

    AttitudeEKF ekf(5 * kArcsec, 30 * kArcsec, DECIMAL(2e-4));
    CHECK_FALSE(ekf.Reset(newer, rate, baseline));
    CHECK_FALSE(ekf.Predict(DECIMAL(1.0)));

    // a quarter turn is the limit, whichever axis it is about
    const Vec3 axes[] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {DECIMAL(0.6), DECIMAL(-0.48), DECIMAL(0.64)}};
    for (const Vec3 &axis : axes) {
        INFO("axis " << axis.x << ", " << axis.y << ", " << axis.z);
        CHECK(ekf.Reset(older, axis * (DECIMAL_M_PI / 2 * DECIMAL(0.99) / baseline), baseline));
        CHECK_FALSE(ekf.Reset(older, axis * (DECIMAL_M_PI / 2 * DECIMAL(1.01) / baseline), baseline));
    }
}

TEST_CASE("EKF covariance is positive definite or the call is refused, for turns about a skew axis", "[kalman]") {
    const Quaternion identity(1, 0, 0, 0);
    // sigmas computed in double and then narrowed, as in the report of this case
    const decimal crossSigma = static_cast<decimal>(5 * (M_PI / 648000.0));
    const decimal rollSigma = static_cast<decimal>(30 * (M_PI / 648000.0));

    // The reported case: in float mode this used to return true with an indefinite covariance.
    {
        AttitudeEKF ekf(crossSigma, rollSigma, DECIMAL(1e-7));
        REQUIRE(ekf.Reset(identity, Vec3(DECIMAL(0.01), DECIMAL(0.02), DECIMAL(0.03)), DECIMAL(0.1)));
        const AttitudeEKF::Mat6 before = ekf.Covariance();
        const Quaternion attitudeBefore = ekf.AttitudeQuaternion();
        if (ekf.Predict(DECIMAL(494.0))) {
            CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
        } else {
            CHECK(ekf.Covariance() == before);
            CHECK(AngleBetween(ekf.AttitudeQuaternion(), attitudeBefore) == 0);
            CHECK((ekf.Rate() - Vec3(DECIMAL(0.01), DECIMAL(0.02), DECIMAL(0.03))).Magnitude() == 0);
        }
    }

    // Around that case: several skew rates, gaps on both sides of where float mode gives out,
    // with and without process noise. Every call either is refused and changes nothing, or
    // leaves a covariance that factors in double precision.
    const Vec3 rates[] = {{DECIMAL(0.01), DECIMAL(0.02), DECIMAL(0.03)}, {DECIMAL(0.02), DECIMAL(-0.01), DECIMAL(0.005)},
                          {DECIMAL(-0.003), DECIMAL(0.004), DECIMAL(0.012)}, {DECIMAL(0.03), DECIMAL(0.001), DECIMAL(-0.001)}};
    int committed = 0, refused = 0;
    for (const Vec3 &rate : rates) {
        for (decimal baseline : {DECIMAL(0.1), DECIMAL(1.0)}) {
            for (decimal accelerationSigma : {DECIMAL(0.0), DECIMAL(1e-7), DECIMAL(1e-5)}) {
                for (int gap = 20; gap <= 1500; gap += 37) {
                    INFO("rate " << rate.x << ", " << rate.y << ", " << rate.z << ", baseline " << baseline
                         << ", acceleration sigma " << accelerationSigma << ", gap " << gap);
                    AttitudeEKF ekf(crossSigma, rollSigma, accelerationSigma);
                    REQUIRE(ekf.Reset(identity, rate, baseline));
                    const AttitudeEKF::Mat6 before = ekf.Covariance();
                    if (!ekf.Predict(static_cast<decimal>(gap))) {
                        refused++;
                        CHECK(ekf.Covariance() == before);
                        continue;
                    }
                    committed++;
                    CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
                    const AttitudeEKF::Mat6 predictedCovariance = ekf.Covariance();
                    if (ekf.Update(ekf.AttitudeQuaternion(), DECIMAL(16.27))) {
                        CHECK(IsPositiveDefinite(ekf.Covariance(), 0));
                    } else {
                        CHECK(ekf.Covariance() == predictedCovariance);
                    }
                }
            }
        }
    }
    CAPTURE(committed, refused);
    CHECK(committed > 200);
#ifdef LOST_FLOAT_MODE
    // float mode gives out somewhere in this range, so both outcomes are exercised
    CHECK(refused > 50);
#endif
}
