#include <cmath>
#include <limits>

#include <catch.hpp>

#include "attitude-utils.hpp"
#include "kalman-filter.hpp"

using namespace lost; // NOLINT

namespace {

const decimal kArcsec = DECIMAL(M_PI) / (180 * 3600);

decimal AngleBetween(const Quaternion &a, const Quaternion &b) {
    return QuaternionToRotationVector(a * b.Conjugate()).Magnitude();
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
    CHECK_FALSE(ekf.Reset(start, rate, std::numeric_limits<decimal>::min()));  // covariance overflows

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

TEST_CASE("EKF with unusable noise settings never accepts a measurement", "[kalman]") {
    const Quaternion identity(1, 0, 0, 0);
    AttitudeEKF zero(0, 0, 0);
    CHECK(zero.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
    CHECK_FALSE(zero.Update(identity, DECIMAL(16.27)));

    AttitudeEKF broken(std::numeric_limits<decimal>::quiet_NaN(), 30 * kArcsec, DECIMAL(0.01));
    CHECK_FALSE(broken.Reset(identity, Vec3(0, 0, 0), DECIMAL(0.1)));
    CHECK_FALSE(broken.Update(identity, DECIMAL(16.27)));
    CHECK(broken.AttitudeQuaternion().IsUnit(kUnitTolerance));
}
