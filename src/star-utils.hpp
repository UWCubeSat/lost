#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <algorithm>
#include <vector>

#include "attitude-utils.hpp"
#include "serialize-helpers.hpp"

namespace lost {

/// A star from the Bright Star Catalog.
class CatalogStar {
public:
    CatalogStar() = default;

    /**
     * Create a CatalogStar using its celestial/spherical coordinates.
     * @param raj2000,dej2000 The right ascension and declination of the star in the year 2000
     * @param magnitude See ::magnitude
     * @param name See ::name
     */
    CatalogStar(decimal raj2000, decimal dej2000, int magnitude, int name) :
        spatial(SphericalToSpatial(raj2000, dej2000)), magnitude(magnitude), name(name) {}

    CatalogStar(Vec3 spatial, int magnitude, int name) :
        spatial(spatial), magnitude(magnitude), name(name) {}

    /**
     * The point on the unit sphere where the star lies.
     * The earth is at the center of the sphere. A star with right ascension and declination both equal to zero will lie at (1,0,0).
     */
    Vec3 spatial;
    /**
     * The magnitude of the star, with the decimal point shifted two places right.
     * I.e., multiply this by 10<sup>-2</sup> to get the magnitude.
     */
    int magnitude;
    /**
     * A unique number which unambiguously identifies the catalog star.
     * The catalog is an array, so in algorithms we usually refer to catalog stars by their index in the catalog. However, if the catalog is filtered (say, to remove stars that are too dim or too close to each other), the indexes change, which makes everything really hard to debug because indexes change based on how the catalog was filtered. That's what names are for: They're fixed and don't change.
     */
    int name;
};

/**
 * A "centroid" detected in an image.
 * This represents a star from an image which has not necessarily been identified.
 */
class Star {
public:
    Star(decimal x, decimal y, decimal radiusX, decimal radiusY, int magnitude) :
        position({x, y}), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) {};

    /// Convenience constructor that sets Star.radiusY = radiusX and Star.magnitude = 0
    Star(decimal x, decimal y, decimal radiusX) : Star(x, y, radiusX, radiusX, 0) {};

    /// Create a zeroed-out star. Fields should be set immediately after construction.
    Star() : Star(DECIMAL(0.0), DECIMAL(0.0), DECIMAL(0.0)) {};

    /// The (x,y) pixel coordinates in the image (top left is 0,0)
    Vec2 position;
    /// Approximate horizontal radius of the bright area in pixels.
    decimal radiusX;
    /// Approximate vertical radius of the bright area in pixels.
    decimal radiusY;
    /**
     * A relative measure of magnitude of the star. Larger is brighter.
     * It's impossible to tell the true magnitude of the star from the image, without really good camera calibration. Anyway, this field is not meant to correspond to the usual measurement of magnitude. Instead, it's just some measure of brightness which may be specific to the centroiding algorithm. For example, it might be the total number of bright pixels in the star.
     */
    int magnitude;
    // eccentricity?
};

/**
 * Records that a certain Star (detected in the image) corresponds to a certain CatalogStar.
 * Only stores indexes into an array of Star objects and an array of CatalogStar objects, which are stored elsewhere.
 */
class StarIdentifier {
public:
    StarIdentifier(int starIndex, int catalogIndex, int weight)
        : starIndex(starIndex), catalogIndex(catalogIndex), weight(weight) { };
    /// Sets StarIdentifier.weight = 1
    StarIdentifier(int starIndex, int catalogIndex)
        : StarIdentifier(starIndex, catalogIndex, DECIMAL(1.0)) { };

    // does not check weight
    bool operator==(const StarIdentifier& other) const {
        return starIndex == other.starIndex &&
            catalogIndex == other.catalogIndex;
    }

    /// An index into an array of Star objects.
    int starIndex;
    /// An index into an array of CatalogStar objects.
    int catalogIndex;
    /// A weight indicating the confidence of this idenification. Often just 1.
    decimal weight;
};

typedef std::vector<CatalogStar> Catalog;
typedef std::vector<Star> Stars;
typedef std::vector<StarIdentifier> StarIdentifiers;

void SerializeCatalog(SerializeContext *, const Catalog &, bool inclMagnitude, bool inclName);
// sets magnited and name to whether the catalog in the database contained magnitude and name
Catalog DeserializeCatalog(DeserializeContext *des, bool *inclMagnitudeReturn, bool *inclNameReturn);
Catalog::const_iterator FindNamedStar(const Catalog &catalog, int name);

/// returns some relative brightness measure, which is proportional to the total number of photons received from a star.
/// As always, the magnitude is actually 100* the usual magnitude
decimal MagToBrightness(int magnitude);

/**
 * Remove unwanted stars from an unfiltered catalog.
 *
 * TODO: Don't necessarily remove both stars when they're within minSeparation. Instead, if the
 * brightnesses are different enough, just keep the brightest one!
 *
 * @param maxMagnitude Should be 100*(magnitude), just like CatalogStar::magnitude. The narrowed
 * catalog will only contain stars at least as bright as that (i.e., lower magnitude). Pass
 * something like 9999 to not filter based on magnitude.
 *
 * @param maxStars No more than maxStars many stars will be in the narrowed catalog (prefers
 * brightest). Pass something like 99999 to not filter maxStars.
 *
 * @param minDistance Any star which is within minDistance of another star is removed from the
 * catalog. Pass something like -1 to not filter based on minDistance. TODO: In the future, some
 * provision may be made so that instead of being removed, these stars get a special distinguishing
 * mark, so that star-id algorithms can reason about them instead of thinking they are false stars,
 * but no such provision is made yet.
 */
Catalog NarrowCatalog(const Catalog &, int maxMagnitude, int maxStars, decimal minSeparation);

// Hash binned edge ratios of a Tetra star pattern into an index in the pattern table.
size_t TetraFindPatternIndex(const int* patternEdgeRatios,
                             size_t numRatios,
                             int numEdgeRatioBins,
                             size_t tableSize);

// Generates all ways to choose k items from n, prioritizing patterns whose
// dimmest star is as bright as possible (intuitively we want to avoid a bad centroid ruining a
// pattern) E.g. for k=4, n=8: 0123, 0124, 0134, 0234, 1234, 0125, ... Formally, all combinations
// using only items 0..m are exhausted before item m+1 appears
class CombinationIterator {
   public:
    CombinationIterator(int k, int n) : k_(k), n_(n), indices_(k), isFirst_(true) {
        for (int i = 0; i < k; i++)
            indices_[i] = i;
    }

    bool next() {
        if (isFirst_) {
            isFirst_ = false;
            return n_ >= k_;
        }
        int i = 0;
        // Find leftmost index that can be incremented
        while (i < k_ && indices_[i] + 1 == (i + 1 < k_ ? indices_[i + 1] : n_)) {
            i++;
        }
        if (i == k_) return false;
        indices_[i]++;
        // Reset everything before the index we just changed
        for (int j = 0; j < i; j++)
            indices_[j] = j;
        return true;
    }

    const std::vector<int>& indices() const { return indices_; }

   private:
    int k_, n_;
    std::vector<int> indices_;
    bool isFirst_;
};

struct TetraConstants {
    static constexpr int numPatternStars = 4;
};

// Check whether all pairs of stars have angular distance within maxFov.
// SpatialAccessor(int i) must return const Vec3&.
template <typename SpatialAccessor>
bool AllPairsWithinFov(int numStars, SpatialAccessor spatial, decimal cosMaxFov) {
    for (int i = 0; i < numStars; i++) {
        for (int j = i + 1; j < numStars; j++) {
            if (spatial(i) * spatial(j) <= cosMaxFov) return false;
        }
    }
    return true;
}

// Compute sorted edge ratios for a star pattern. SpatialAccessor(int i) must return const Vec3&.
template <typename SpatialAccessor>
std::vector<decimal> TetraEdgeRatios(int numStars, SpatialAccessor spatial) {
    std::vector<decimal> edgeLengths;
    for (int i = 0; i < numStars; i++) {
        for (int j = i + 1; j < numStars; j++) {
            edgeLengths.push_back((spatial(j) - spatial(i)).Magnitude());
        }
    }
    std::sort(edgeLengths.begin(), edgeLengths.end());
    decimal longestEdgeLen = edgeLengths.back();
    std::vector<decimal> edgeRatios;
    for (size_t i = 0; i < edgeLengths.size() - 1; i++) {
        edgeRatios.push_back(edgeLengths[i] / longestEdgeLen);
    }
    return edgeRatios;
}

}  // namespace lost

#endif
