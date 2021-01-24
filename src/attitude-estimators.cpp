#include "attitude-estimators.hpp"

namespace lost {
    Quaternion AttitudeEstimationAlgorithm::Go(const Camera &, const Stars &, const Catalog &, const StarIdentifiers &) {
        return Quaternion();
    }
}