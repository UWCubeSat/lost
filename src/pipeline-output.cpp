#include "pipeline-output.hpp"

#include <cairo/cairo.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>

#include "attitude-utils.hpp"
#include "cairo-utils.hpp"
#include "decimal.hpp"
#include "io-util.hpp"
#include "pipeline-input.hpp"
#include "star-utils.hpp"

namespace lost {

////////////////
// COMPARISON //
////////////////

/**
 * The result of comparing actual and expected centroids
 * Used for debugging and benchmarking
 */
class CentroidComparison {
public:
    CentroidComparison() : meanError(0.0), numCorrectCentroids(0), numExtraCentroids(0) { }
    /**
     * Average distance from actual to expected centroids (in pixels)
     * Only correct centroids are considered in this average.
     */
    decimal meanError;

    /**
     * Number of actual stars within the centroiding threshold of an expected star.
     */
    decimal numCorrectCentroids;

    /**
     * Stars in actual but not expected. Ideally 0
     * This is a decimal because we may average multiple centroid comparisons together.
     */
    decimal numExtraCentroids;

    // We no longer have a num missing because often generated stars have too low of a signal-to-noise ratio and the centroid algo won't pick them up.
};

/// Create a mapping, where keys are indices into `one` and values are indices of all centroids in
/// `two` whose distance to the current `one` star is <=threshold. In a (ordered) multimap, the
/// insertion order is preserved for elements with the same key, and indeed we'll sort the elements
/// corresponding to each key by distance from the corresponding `one` star.
static std::multimap<int, int> FindClosestCentroids(decimal threshold,
                                                    const Stars &one,
                                                    const Stars &two) {
    std::multimap<int, int> result;

    for (int i = 0; i < (int)one.size(); i++) {
        std::vector<std::pair<decimal, int>> closest;
        for (int k = 0; k < (int)two.size(); k++) {
            decimal currDistance = (one[i].position - two[k].position).Magnitude();
            if (currDistance <= threshold) {
                closest.emplace_back(currDistance, k);
            }
        }
        std::sort(closest.begin(), closest.end());
        for (const std::pair<decimal, int> &pair : closest) {
            result.emplace(i, pair.second);
        }
    }

    return result;
}

/**
 * Compare expected and actual centroids.
 * Useful for debugging and benchmarking.
 * @param threshold The maximum number of pixels apart two centroids can be to be considered the same.
 */
static CentroidComparison CentroidsCompare(decimal threshold,
                                    const Stars &expected,
                                    const Stars &actual) {

    // TODO: Somehow penalize when multiple centroids correspond to the same expected star (i.e.,
    // one star turned into multiple centroids). That should probably be considered an extra
    // centroid, but rn it isn't.

    CentroidComparison result;
    // maps from indexes in each list to the closest centroid from other list
    std::multimap<int, int> actualToExpected = FindClosestCentroids(threshold, actual, expected);

    for (int i = 0; i < (int)actualToExpected.size(); i++) {
        auto closest = actualToExpected.find(i);
        if (closest == actualToExpected.end()) {
            result.numExtraCentroids++;
        } else {
            result.meanError += (actual[i].position - expected[closest->second].position).Magnitude();
            result.numCorrectCentroids++;
        }
    }
    result.meanError /= result.numCorrectCentroids;

    return result;
}

static CentroidComparison CentroidComparisonsCombine(std::vector<CentroidComparison> comparisons) {
    assert(comparisons.size() > 0);

    CentroidComparison result;

    for (const CentroidComparison &comparison : comparisons) {
        result.meanError += comparison.meanError;
        result.numCorrectCentroids += comparison.numCorrectCentroids;
        result.numExtraCentroids += comparison.numExtraCentroids;
    }

    result.meanError /= comparisons.size();
    result.numExtraCentroids /= comparisons.size();
    result.numCorrectCentroids /= comparisons.size();

    return result;
}

// (documentation in hpp)
StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                // use these to map indices to names for the respective lists of StarIdentifiers
                                const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                decimal centroidThreshold,
                                const Stars &expectedStars, const Stars &inputStars) {

    StarIdComparison result = {
        0, // correct
        0, // incorrect
        0, // total
    };

    // EXPECTED STAR IDS

    // map from expected star indices to expected catalog indices (basically flattening the expected star-ids)
    std::vector<int> expectedCatalogIndices(expectedStars.size(), -1);
    for (const StarIdentifier &starId : expected) {
        assert(0 <= starId.starIndex && starId.starIndex <= (int)expectedStars.size());
        assert(0 <= starId.catalogIndex && starId.catalogIndex <= (int)expectedCatalog.size());
        expectedCatalogIndices[starId.starIndex] = starId.catalogIndex;
    }

    // FIND NEAREST CENTROIDS

    std::multimap<int, int> inputToExpectedCentroids = FindClosestCentroids(centroidThreshold, inputStars, expectedStars);
    // std::multimap<int, int> expectedToInputCentroids = FindClosestCentroids(centroidThreshold, expectedStars, inputStars);

    // COMPUTE TOTAL
    // Count the number of expected stars with at least one input star near them
    for (int i = 0; i < (int)inputStars.size(); i++) {
        // make sure there's at least one expected star near this input star which has an identification
        auto closestRange = inputToExpectedCentroids.equal_range(i);
        bool found = false;
        for (auto it = closestRange.first; it != closestRange.second; it++) {
            if (expectedCatalogIndices[it->second] != -1) {
                found = true;
                break;
            }
        }
        if (found) {
            result.numTotal++;
        }
    }

    // COMPUTE CORRECT AND INCORRECT

    std::vector<bool> identifiedInputCentroids(inputStars.size(), false);
    for (const StarIdentifier &starId : actual) {
        // as later, there shouldn't be duplicate starIndex. This indicates a bug in the star-id algorithm, not comparison code.
        assert(!identifiedInputCentroids[starId.starIndex]);
        identifiedInputCentroids[starId.starIndex] = true;
        assert(0 <= starId.starIndex && starId.starIndex <= (int)inputStars.size());
        assert(0 <= starId.catalogIndex && starId.catalogIndex <= (int)actualCatalog.size());

        // Check that there's at least one expected centroid in range which agrees with your identification.
        auto expectedCentroidsInRange = inputToExpectedCentroids.equal_range(starId.starIndex);
        bool found = false;
        for (auto it = expectedCentroidsInRange.first; it != expectedCentroidsInRange.second; it++) {
            int expectedCatalogIndex = expectedCatalogIndices[it->second];
            if (expectedCatalogIndex != -1
                && expectedCatalog[expectedCatalogIndex].name == actualCatalog[starId.catalogIndex].name) {

                result.numCorrect++;
                found = true;
                break;
            }
        }

        // Either there's no expected centroid in range, or none of them agree with the identification.
        if (!found) {
            result.numIncorrect++;
        }
    }

    return result;
}

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

typedef void (*PipelineComparator)(std::ostream &os,
                                   const PipelineInputList &,
                                   const std::vector<PipelineOutput> &,
                                   const PipelineOptions &);

/// Plotter suitable for `cairo_surface_write_to_png_stream` which simply writes to an std::ostream
static cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

/// Plots the input image with no annotation to `os`
static void PipelineComparatorPlotRawInput(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &,
                                    const PipelineOptions &) {

    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Plot the annotated input image to `os`
// TODO: should probably use Expected methods, not Input methods, because future PipelineInputs could add noise to the result of the Input methods.
static void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &,
                                 const PipelineOptions &) {
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    assert(expected[0]->InputStars() != NULL);
    SurfacePlot("pipeline input",
                cairoSurface,
                *expected[0]->InputStars(),
                expected[0]->InputStarIds(),
                &expected[0]->GetCatalog(),
                expected[0]->InputAttitude(),
                // green
                0.0, 1.0, 0.0, 0.6);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

static void PipelineComparatorPlotExpected(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &,
                                    const PipelineOptions &) {
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    assert(expected[0]->ExpectedStars() != NULL);
    SurfacePlot("expected output",
                cairoSurface,
                *expected[0]->ExpectedStars(),
                expected[0]->ExpectedStarIds(),
                &expected[0]->GetCatalog(),
                expected[0]->ExpectedAttitude(),
                // blu
                0.2, 0.5, 1.0, 0.7);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Compare the actual and expected centroids, printing key stats to `os`
static void PipelineComparatorCentroids(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual,
                                 const PipelineOptions &values) {
    int size = (int)expected.size();

    decimal threshold = values.centroidCompareThreshold;

    std::vector<CentroidComparison> comparisons;
    for (int i = 0; i < size; i++) {
        comparisons.push_back(CentroidsCompare(threshold,
                                               *(expected[i]->ExpectedStars()),
                                               *(actual[i].stars)));
    }

    CentroidComparison result = CentroidComparisonsCombine(comparisons);
    os << "centroids_num_correct " << result.numCorrectCentroids << std::endl
       << "centroids_num_extra " << result.numExtraCentroids << std::endl
       << "centroids_mean_error " << result.meanError << std::endl;
}

static void PrintCentroids(const std::string &prefix,
                           std::ostream &os,
                           const Catalog &catalog,
                           const std::vector<Stars> &starses,
                           // May be NULL. Should be the only the first starId, because we don't have any reasonable aggregative action to perform.
                           const StarIdentifiers *starIds) {
    assert(starses.size() > 0);
    decimal avgNumStars = 0;
    for (const Stars &stars : starses) {
        avgNumStars += stars.size();
    }
    avgNumStars /= starses.size();

    os << "num_" << prefix << "_centroids " << avgNumStars << std::endl;
    if (starses.size() == 1) {
        const Stars &stars = starses[0];
        for (int i = 0; i < (int)stars.size(); i++) {
            os << prefix << "_centroid_" << i << "_x " << stars[i].position.x << std::endl;
            os << prefix << "_centroid_" << i << "_y " << stars[i].position.y << std::endl;
            if (starIds) {
                for (const StarIdentifier &starId : *starIds) {
                    if (starId.starIndex == i) {
                        os << prefix << "_centroid_" << i << "_id " << catalog[starId.catalogIndex].name << std::endl;
                    }
                }
            }
        }
    }
}

/// Print a list of centroids to `os`
static void PipelineComparatorPrintExpectedCentroids(std::ostream &os,
                                                     const PipelineInputList &expected,
                                                     const std::vector<PipelineOutput> &, // actual
                                                     const PipelineOptions &) {
    assert(expected.size() > 0);
    assert(expected[0]->ExpectedStars());

    std::vector<Stars> expectedStarses;
    for (const auto &input : expected) {
        expectedStarses.push_back(*input->ExpectedStars());
    }
    PrintCentroids("expected",
                   os,
                   expected[0]->GetCatalog(),
                   expectedStarses,
                   expected[0]->ExpectedStarIds());
}

static void PipelineComparatorPrintInputCentroids(std::ostream &os,
                                                  const PipelineInputList &expected,
                                                  const std::vector<PipelineOutput> &, // actual
                                                  const PipelineOptions &) {
    assert(expected.size() > 0);
    assert(expected[0]->InputStars());

    std::vector<Stars> inputStarses;
    for (const auto &input : expected) {
        inputStarses.push_back(*input->InputStars());
    }
    PrintCentroids("input",
                   os,
                   expected[0]->GetCatalog(),
                   inputStarses,
                   expected[0]->InputStarIds());
}

static void PipelineComparatorPrintActualCentroids(std::ostream &os,
                                                   const PipelineInputList &expected, // expected
                                                   const std::vector<PipelineOutput> &actual,
                                                   const PipelineOptions &values) {
    assert(actual.size() > 0);
    assert(actual[0].stars);

    std::vector<Stars> actualStarses;
    for (const auto &output : actual) {
        actualStarses.push_back(*output.stars);
    }
    PrintCentroids("actual",
                   os,
                   actual[0].catalog,
                   actualStarses,
                   actual[0].starIds.get());

    if (expected.size() == 1 && expected[0]->ExpectedStars() && expected[0]->ExpectedStarIds()) {
        // also print expected ID of each one
        const Stars &actualStars = *actual[0].stars;
        const Stars &expectedStars = *expected[0]->ExpectedStars();
        std::multimap<int, int> actualToExpectedCentroids = FindClosestCentroids(values.centroidCompareThreshold, actualStars, expectedStars);
        for (int i = 0; i < (int)actualStars.size(); i++) {
            auto range = actualToExpectedCentroids.equal_range(i);
            auto it = range.first;
            auto end = range.second;
            for (int j = 0;
                 it != end;
                 j++, it++) {

                int expectedCentroidIndex = it->second;
                bool foundIt = false; // just to be sure
                for (const StarIdentifier &starId : *expected[0]->ExpectedStarIds()) {
                    if (starId.starIndex == expectedCentroidIndex) {
                        assert(!foundIt);
                        int expectedName = expected[0]->GetCatalog()[starId.catalogIndex].name;
                        std::cout << "actual_centroid_" << i << "_expected_id_" << j << " " << expectedName << std::endl;
                        foundIt = true;
                    }
                }
            }
        }
    }
}

/// Plot an annotated image where centroids are annotated with their centroid index. For debugging.
/// Use whatever stars were input into the star-id algo (so either actual centroids, or inputstars)
static void PipelineComparatorPlotCentroidIndices(std::ostream &os,
                                           const PipelineInputList &expected,
                                           const std::vector<PipelineOutput> &actual,
                                           const PipelineOptions &) {
    const Stars &stars = actual[0].stars ? *actual[0].stars : *expected[0]->InputStars();
    StarIdentifiers identifiers;
    for (int i = 0; i < (int)stars.size(); i++) {
        identifiers.push_back(StarIdentifier(i, i));
    }
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot("centroid indices (input)",
                cairoSurface,
                stars,
                &identifiers,
                &actual[0].catalog,
                NULL,
                // orange
                1.0, 0.5, 0.0, 0.5,
                // don't resolve names
                true);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Plot the image annotated with output data computed by the star tracking algorithms.
static void PipelineComparatorPlotOutput(std::ostream &os,
                                         const PipelineInputList &expected,
                                         const std::vector<PipelineOutput> &actual,
                                         const PipelineOptions &) {
    // don't need to worry about mutating the surface; InputImageSurface returns a fresh one
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot("pipeline output",
                cairoSurface,
                actual[0].stars ? *actual[0].stars : *expected[0]->InputStars(),
                actual[0].starIds.get(),
                &actual[0].catalog,
                actual[0].attitude.get(),
                // red
                1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Compare the expected and actual star identifiers.
static void PipelineComparatorStarIds(std::ostream &os,
                                      const PipelineInputList &expected,
                                      const std::vector<PipelineOutput> &actual,
                                      const PipelineOptions &values) {
    int numImagesCorrect = 0;
    int numImagesIncorrect = 0;
    int numImagesTotal = expected.size();
    for (int i = 0; i < numImagesTotal; i++) {
        // since the actual star IDs exist, it must have gotten input from somewhere!
        // TODO: overhaul: It seems that in these comparators there should be a more fundamental way to figure out the input that was actually sent to a stage.
        // I.e., instead of having expected and actual arguments, have some sort of PipelineRunSummary object, where the InputStars method looks at actual, then input.
        assert(actual[i].stars.get() || expected[i]->InputStars());

        const Stars &inputStars = actual[i].stars.get()
            ? *actual[i].stars.get()
            : *expected[i]->InputStars();
        StarIdComparison comparison =
            StarIdsCompare(*expected[i]->ExpectedStarIds(), *actual[i].starIds,
                           expected[i]->GetCatalog(), actual[i].catalog,
                           values.centroidCompareThreshold, *expected[i]->ExpectedStars(), inputStars);

        if (numImagesTotal == 1) {
            os << "starid_num_correct " << comparison.numCorrect << std::endl;
            os << "starid_num_incorrect " << comparison.numIncorrect << std::endl;
            os << "starid_num_total " << comparison.numTotal << std::endl;
        }

        if (comparison.numCorrect > 0 && comparison.numIncorrect == 0) {
            numImagesCorrect++;
        }
        if (comparison.numIncorrect > 0) {
            numImagesIncorrect++;
        }
    }

    // A "correct" image is one where at least two stars are correctly id'd and none are incorrectly id'd
    os << "starid_num_images_correct " << numImagesCorrect << std::endl;
    os << "starid_num_images_incorrect " << numImagesIncorrect << std::endl;
}

static void PrintAttitude(std::ostream &os, const std::string &prefix, const Attitude &attitude) {
    if (attitude.IsKnown()) {
        os << prefix << "attitude_known 1" << std::endl;

        EulerAngles spherical = attitude.ToSpherical();
        os << prefix << "attitude_ra " << RadToDeg(spherical.ra) << std::endl;
        os << prefix << "attitude_de " << RadToDeg(spherical.de) << std::endl;
        os << prefix << "attitude_roll " << RadToDeg(spherical.roll) << std::endl;

        Quaternion q = attitude.GetQuaternion();
        os << prefix << "attitude_i " << q.i << std::endl;
        os << prefix << "attitude_j " << q.j << std::endl;
        os << prefix << "attitude_k " << q.k << std::endl;
        os << prefix << "attitude_real " << q.real << std::endl;

    } else {
        os << prefix << "attitude_known 0" << std::endl;
    }
}

/// Print the identifed attitude to `os` in Euler angle format.
static void PipelineComparatorPrintAttitude(std::ostream &os,
                                            const PipelineInputList &,
                                            const std::vector<PipelineOutput> &actual,
                                            const PipelineOptions &) {
    assert(actual.size() == 1);
    assert(actual[0].attitude);
    PrintAttitude(os, "", *actual[0].attitude);
}

static void PipelineComparatorPrintExpectedAttitude(std::ostream &os,
                                                   const PipelineInputList &expected,
                                                   const std::vector<PipelineOutput> &,
                                                   const PipelineOptions &) {
    assert(expected.size() == 1);
    assert(expected[0]->ExpectedAttitude());
    PrintAttitude(os, "expected_", *expected[0]->ExpectedAttitude());
}

/// Compare the actual and expected attitudes.
static void PipelineComparatorAttitude(std::ostream &os,
                                       const PipelineInputList &expected,
                                       const std::vector<PipelineOutput> &actual,
                                       const PipelineOptions &values) {

    // TODO: use Wahba loss function (maybe average per star) instead of just angle. Also break
    // apart roll error from boresight error. This is just quick and dirty for testing

    decimal angleThreshold = DegToRad(values.attitudeCompareThreshold);

    decimal attitudeErrorSum = 0.0f;
    int numCorrect = 0;
    int numIncorrect = 0;

    for (int i = 0; i < (int)expected.size(); i++) {
        if (actual[i].attitude->IsKnown()) {
            Quaternion expectedQuaternion = expected[i]->ExpectedAttitude()->GetQuaternion();
            Quaternion actualQuaternion = actual[i].attitude->GetQuaternion();
            decimal attitudeError = (expectedQuaternion * actualQuaternion.Conjugate()).SmallestAngle();
            assert(attitudeError >= 0);

            if (attitudeError <= angleThreshold) {
                attitudeErrorSum += attitudeError;
                numCorrect++;
            } else {
                numIncorrect++;
            }
        }
    }

    decimal attitudeErrorMean = DECIMAL(attitudeErrorSum) / numCorrect;
    decimal fractionCorrect = DECIMAL(numCorrect) / expected.size();
    decimal fractionIncorrect = DECIMAL(numIncorrect) / expected.size();

    os << "attitude_error_mean " << attitudeErrorMean << std::endl;
    os << "attitude_availability " << fractionCorrect << std::endl;
    os << "attitude_error_rate " << fractionIncorrect << std::endl;
}

static void PrintTimeStats(std::ostream &os, const std::string &prefix, const std::vector<long long> &times) {
    assert(times.size() > 0);

    // print average, min, max, and 95% max
    long long sum = 0;
    long long min = LONG_MAX;
    long long max = 0;
    for (int i = 0; i < (int)times.size(); i++) {
        assert(times[i] > 0);
        sum += times[i];
        min = std::min(min, times[i]);
        max = std::max(max, times[i]);
    }
    long average = sum / times.size();
    std::vector<long long> sortedTimes = times;
    std::sort(sortedTimes.begin(), sortedTimes.end());
    // what really is the 95th percentile? Being conservative, we want to pick a value that at least
    // 95% of the times are less than. This means: (1) finding the number of times, (2) Finding
    // Math.ceil(0.95 * numTimes), and (3) subtracting 1 to get the index.
    int ninetyFiveIndex = (int)std::ceil(0.95 * times.size()) - 1;
    assert(ninetyFiveIndex >= 0);
    long long ninetyFifthPercentile = sortedTimes[ninetyFiveIndex];

    os << prefix << "_average_ns " << average << std::endl;
    os << prefix << "_min_ns " << min << std::endl;
    os << prefix << "_max_ns " << max << std::endl;
    os << prefix << "_95%_ns " << ninetyFifthPercentile << std::endl;
}

/// For each stage of the pipeline, print statistics about how long it took to run.
static void PipelineComparatorPrintSpeed(std::ostream &os,
                                    const PipelineInputList &,
                                    const std::vector<PipelineOutput> &actual,
                                    const PipelineOptions &) {
    std::vector<long long> centroidingTimes;
    std::vector<long long> starIdTimes;
    std::vector<long long> attitudeTimes;
    std::vector<long long> totalTimes;
    for (int i = 0; i < (int)actual.size(); i++) {
        long long totalTime = 0;
        if (actual[i].centroidingTimeNs > 0) {
            centroidingTimes.push_back(actual[i].centroidingTimeNs);
            totalTime += actual[i].centroidingTimeNs;
        }
        if (actual[i].starIdTimeNs > 0) {
            starIdTimes.push_back(actual[i].starIdTimeNs);
            totalTime += actual[i].starIdTimeNs;
        }
        if (actual[i].attitudeEstimationTimeNs > 0) {
            attitudeTimes.push_back(actual[i].attitudeEstimationTimeNs);
            totalTime += actual[i].attitudeEstimationTimeNs;
        }
        totalTimes.push_back(totalTime);
    }
    if (centroidingTimes.size() > 0) {
        PrintTimeStats(os, "centroiding", centroidingTimes);
    }
    if (starIdTimes.size() > 0) {
        PrintTimeStats(os, "starid", starIdTimes);
    }
    if (attitudeTimes.size() > 0) {
        PrintTimeStats(os, "attitude", attitudeTimes);
    }
    if (centroidingTimes.size() > 0 || starIdTimes.size() > 0 || attitudeTimes.size() > 0) {
        PrintTimeStats(os, "total", totalTimes);
    }
}

/**
 * Print or otherwise analyze the results of (perhaps multiple) runs of a star tracking pipeline.
 * Uses the command line options in `values` to determine which analyses to run. Examples include plotting an annotated output image to a png file, comparing the actual and expected centroids, etc
 */
void PipelineComparison(const PipelineInputList &expected,
                        const std::vector<PipelineOutput> &actual,
                        const PipelineOptions &values) {
    if (actual.size() == 0) {
        std::cerr << "ERROR: No output! Did you specify any input images? Try --png or --generate." << std::endl;
        exit(1);
    }

    assert(expected.size() == actual.size() && expected.size() > 0);

    // TODO: Remove the asserts and print out more reasonable error messages.

#define LOST_PIPELINE_COMPARE(precondition, errmsg, comparator, path, isBinary) do { \
        if (precondition) {                                             \
            UserSpecifiedOutputStream pos(path, isBinary);              \
            comparator(pos.Stream(), expected, actual, values);         \
        } else {                                                        \
            std::cerr << "ERROR: Comparator not applicable: " << errmsg << std::endl; \
            exit(1);                                                    \
        }                                                               \
    } while (0)

    if (values.plotRawInput != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1,
                              "--plot-raw-input requires exactly 1 input image, but " + std::to_string(expected.size()) + " many were provided.",
                              PipelineComparatorPlotRawInput, values.plotRawInput, true);
    }

    if (values.plotInput != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1 && expected[0]->InputStars(),
                              "--plot-input requires exactly 1 input image, and for centroids to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPlotInput, values.plotInput, true);
    }
    if (values.plotExpected != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1 && expected[0]->ExpectedStars(),
                              "--plot-expected-input requires exactly 1 input image, and for expected centroids to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPlotExpected, values.plotExpected, true);
    }
    if (values.plotOutput != "") {
        LOST_PIPELINE_COMPARE(actual.size() == 1 && (actual[0].stars || actual[0].starIds),
                              "--plot-output requires exactly 1 output image, and for either centroids or star IDs to be available on that output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPlotOutput, values.plotOutput, true);
    }
    if (values.printExpectedCentroids != "") {
        LOST_PIPELINE_COMPARE(expected[0]->ExpectedStars(),
                              "--print-expected-centroids requires at least 1 input with expected centroids. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPrintExpectedCentroids, values.printExpectedCentroids, false);
    }
    if (values.printInputCentroids != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputStars(),
                              "--print-input-centroids requires at least 1 input with centroids. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPrintInputCentroids, values.printInputCentroids, false);
    }
    if (values.printActualCentroids != "") {
        LOST_PIPELINE_COMPARE(actual[0].stars,
                              "--print-actual-centroids requires at least 1 output image, and for centroids to be available on the output images. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintActualCentroids, values.printActualCentroids, false);
    }
    if (values.plotCentroidIndices != "") {
        LOST_PIPELINE_COMPARE(expected.size() == 1 && expected[0]->InputImage(),
                              "--plot-centroid-indices requires exactly 1 input with image. " + std::to_string(expected.size()) + " many inputs were provided.",
                              PipelineComparatorPlotCentroidIndices, values.plotCentroidIndices, true);
    }
    if (values.compareCentroids != "") {
        LOST_PIPELINE_COMPARE(actual[0].stars && expected[0]->ExpectedStars() && values.centroidCompareThreshold,
                              "--compare-centroids requires at least 1 output image, and for expected centroids to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorCentroids, values.compareCentroids, false);
    }
    if (values.compareStarIds != "") {
        LOST_PIPELINE_COMPARE(expected[0]->ExpectedStarIds() && actual[0].starIds && expected[0]->ExpectedStars(),
                              "--compare-star-ids requires at least 1 output image, and for expected star IDs and centroids to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorStarIds, values.compareStarIds, false);
    }
    if (values.printAttitude != "") {
        LOST_PIPELINE_COMPARE(actual[0].attitude && actual.size() == 1,
                              "--print-attitude requires exactly 1 output image, and for attitude to be available on that output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintAttitude, values.printAttitude, false);
    }
    if (values.printExpectedAttitude != "") {
        LOST_PIPELINE_COMPARE(expected[0]->ExpectedAttitude() && expected.size() == 1,
                              "--print-expected-attitude requires exactly 1 input image, and for expected attitude to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPrintExpectedAttitude, values.printExpectedAttitude, false);
    }
    if (values.compareAttitudes != "") {
        LOST_PIPELINE_COMPARE(actual[0].attitude && expected[0]->ExpectedAttitude() && values.attitudeCompareThreshold,
                              "--compare-attitudes requires at least 1 output image, and for expected attitude to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorAttitude, values.compareAttitudes, false);
    }
    if (values.printSpeed != "") {
        LOST_PIPELINE_COMPARE(actual.size() > 0,
                              // I don't think this should ever actually happen??
                              "--print-speed requires at least 1 output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintSpeed, values.printSpeed, false);
    }

#undef LOST_PIPELINE_COMPARE
}

} // namespace lost
