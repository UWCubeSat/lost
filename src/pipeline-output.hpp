#ifndef PIPELINE_OUTPUT_H
#define PIPELINE_OUTPUT_H

#include <vector>
#include <iostream>

#include "star-utils.hpp"
#include "attitude-utils.hpp"
#include "pipeline-input.hpp"

namespace lost {

/// The result of comparing an actual star identification with the true star idenification, used for testing and benchmarking.
struct StarIdComparison {
    /// The number of centroids in the image which are close to an expected centroid that had an
    /// expected identification the same as the actual identification.
    int numCorrect;

    /// The number of centroids which were either:
    /// + False, but identified as something anyway.
    /// + True, with an identification that did not agree with any sufficiently close expected centroid's expected identification.
    int numIncorrect;

    /// The number of centroids sufficiently close to a true expected star.
    int numTotal;
};

// TODO: rename. Do something with the output
void PipelineComparison(const PipelineInputList &expected,
                        const std::vector<PipelineOutput> &actual,
                        const PipelineOptions &values);

/**
 * Compare expected and actual star identifications.
 * Useful for debugging and benchmarking.
 *
 * The following description is compatible with, but more actionable than, the definitions in the
 * documentation for StarIdComparison. A star-id is *correct* if the centroid is the closest
 * centroid to some expected centroid, and the referenced catalog star is the same one as in the
 * expected star-ids for that centroid. Also permissible is if the centroid is not the closest to
 * any expected centroid, but it has the same star-id as another star closer to the closest expected
 * centroid. All other star-ids are *incorrect* (because they are either identifying false stars, or
 * are incorrect identifications on true stars)
 *
 * The "total" in the result is just the number of input stars.
 */
StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                // use these to map indices to names for the respective lists of StarIdentifiers
                                const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                decimal centroidThreshold,
                                const Stars &expectedStars, const Stars &inputStars);

}

#endif
