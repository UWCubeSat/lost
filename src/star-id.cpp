#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>

#include "star-id.hpp"
#include "databases.hpp"
#include "attitude-utils.hpp"

namespace lost {

StarIdentifiers DummyStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers result;

    for (int i = 0; i < (int)stars.size(); i++) {
        if (rand() > RAND_MAX/5) {
            result.push_back(StarIdentifier(i, rand() % 2000));
        }
    }

    return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    for (int i = 0; i < (int)stars.size(); i++) {  
        std::vector<int16_t> votes(catalog.size(), 0);
        Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                // TODO: find a faster way to do this:
                std::vector<bool> votedInPair(catalog.size(), false);
                Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
                float greatCircleDistance = AngleUnit(iSpatial, jSpatial);
                //give a greater range for min-max Query for bigger radius (GreatCircleDistance)
                float lowerBoundRange = greatCircleDistance - tolerance;
                float upperBoundRange = greatCircleDistance + tolerance;
                long numReturnedPairs;
                const int16_t *lowerBoundSearch = vectorDatabase.FindPairsLiberal(
                    lowerBoundRange, upperBoundRange, &numReturnedPairs);
                //loop from lowerBoundSearch till numReturnedPairs, add one vote to each star in the pairs in the datastructure
                for (const int16_t *k = lowerBoundSearch; k < lowerBoundSearch + numReturnedPairs * 2; k++) {
                    if ((k - lowerBoundSearch) % 2 == 0) {
                        float actualAngle = AngleUnit(catalog[*k].spatial, catalog[*(k+1)].spatial);
                        assert(actualAngle <= greatCircleDistance + tolerance * 2);
                        assert(actualAngle >= greatCircleDistance - tolerance * 2);
                    }
                    if (!votedInPair[*k] || true) {
                        // if (i == 542 && *k == 9085) {
                        //     printf("INC, distance %f from query %f to %f\n", greatCircleDistance,
                        //         lowerBoundRange, upperBoundRange);
                        // }
                        votes[*k]++;
                        votedInPair[*k] = true;
                    }
                }
                // US voting system
            }
        }
        // Find star w most votes
        int16_t maxVotes = votes[0];
        int indexOfMax = 0;
        for (int v = 1; v < (int)votes.size(); v++) {
            if (votes[v] > maxVotes) {
                maxVotes = votes[v];
                indexOfMax = v;
            }
        }
        // if (i == 542) {
        //     for (float dist : vectorDatabase.StarDistances(9085, catalog)) {
        //         printf("Actual 9085 distance: %f\n", dist);
        //     }
        //     puts("Debug star.");
        //     for (int i = 0; i < (int)votes.size(); i++) {
        //         if (votes[i] > maxVotes/2) {
        //             printf("Star %4d received %d votes.\n", catalog[i].name, votes[i]);
        //         }
        //     }
        //     printf("Debug star: Actually voted for %d with %d votes\n",
        //            catalog[indexOfMax].name, maxVotes);
        // }
        // printf("Max votes: %d\n", maxVotes);
        //starIndex = i, catalog index = indexOfMax
        StarIdentifier newStar(i, indexOfMax);
        // Set identified[i] to value of catalog index of star w most votesr
        identified.push_back(newStar);
    }
    //optimizations? N^2
    //https://www.researchgate.net/publication/3007679_Geometric_voting_algorithm_for_star_trackers
    //
    // Do we have a metric for localization uncertainty? Star brighntess?
    //loop i from 1 through n
    std::vector<int16_t> verificationVotes(identified.size(), 0);
    for (int i = 0; i < (int)identified.size(); i++) {
        //loop j from i+1 through n 
        for (int j = i + 1; j < (int)identified.size(); j++) {
            // Calculate distance between catalog stars
            CatalogStar first = catalog[identified[i].catalogIndex];
            CatalogStar second = catalog[identified[j].catalogIndex];
            float cDist = AngleUnit(first.spatial, second.spatial);

            Star firstIdentified = stars[identified[i].starIndex];
            Star secondIdentified = stars[identified[j].starIndex];
            Vec3 firstSpatial = camera.CameraToSpatial(firstIdentified.position);
            Vec3 secondSpatial = camera.CameraToSpatial(secondIdentified.position);
            float sDist = Angle(firstSpatial, secondSpatial);
            
            //if sDist is in the range of (distance between stars in the image +- R)
            //add a vote for the match
            if (abs(sDist - cDist) < tolerance) {
                verificationVotes[i]++;
                verificationVotes[j]++;
            }
        }
    }
    // Find star w most votes
    int maxVotes = verificationVotes.size() > 0 ? verificationVotes[0] : 0;
    for (int v = 1; v < (int)verificationVotes.size(); v++) {
        if (verificationVotes[v] > maxVotes) {
            maxVotes = verificationVotes[v];
        }
    }

    // If the stars are within a certain range of the maximal number of votes, 
    // we consider it correct.
    // maximal votes = maxVotes
    StarIdentifiers verified;
    int thresholdVotes = maxVotes * 3 / 4;
    printf("Verification threshold: %d\n", thresholdVotes);
    for (int i = 0; i < (int)verificationVotes.size(); i++) {
        if (verificationVotes[i] > thresholdVotes) {
            verified.push_back(identified[i]);
        }
    }

    return verified;
}

    /**
     * Strategies:
     * 
     * 1. For each star, enumerate all stars which have the same combination of distances to some
     *  other stars, getting down to a hopefully small (<10) list of candidates for each star, then
     *  do a quad-nested loop to correlate them.
     *
     * 2. Loop through all possible stars in the catalog for star i. Then look at edge ij, using
     * this to select possible j-th stars. If ever there is not a possible j-th star, continue the
     * i-loop. When a possible ij combination is found, loop through k stars according to ik. IF
     * none are found, continue the outer i loop. If some are found, check jk for each one. For each possible ijk triangle, 
     */

class PairDistanceInvolvingIterator {
public:
    // unqualified constructor makes a "past-the-end" iterator
    PairDistanceInvolvingIterator()
        : pairs(NULL), pastTheEnd(NULL) { };

    PairDistanceInvolvingIterator(const int16_t *pairs, long numPairs, int16_t involving)
        : pairs(pairs), pastTheEnd(pairs + numPairs*2), involving(involving) {

        forwardUntilInvolving();
    };

    // PairDistanceInvolvingIterator operator++() {
    //     PairDistanceInvolvingIterator result(*this);
    //     ++(*this);
    //     return result;
    // }

    PairDistanceInvolvingIterator &operator++() {
        assert(hasValue());
        pairs += 2;
        forwardUntilInvolving();
        return *this;
    }

    int16_t operator*() const {
        return curValue;
    }

    bool hasValue() {
        return pairs != pastTheEnd;
    }

    // bool operator==(const PairDistanceInvolvingIterator &other) const {
    //     return ()other.pairs == pairs;
    // }

    // bool operator!=(const PairDistanceInvolvingIterator &other) const {
    //     return !(*this == other);
    // }
private:
    const int16_t *pairs;
    const int16_t *pastTheEnd;
    int16_t involving;
    int16_t curValue;

    // like postfix++, except it's a no-op if already on a valid spot.
    void forwardUntilInvolving() {
        while (pairs != pastTheEnd) {
            if (pairs[0] == involving) {
                curValue = pairs[1];
                return;
            }
            if (pairs[1] == involving) {
                curValue = pairs[0];
                return;
            }
            pairs += 2;
        }
    }
};

void PyramidIdentifyRemainingStars(StarIdentifiers *identifiers,
                                   const Stars &stars,
                                   const Catalog &catalog,
                                   const PairDistanceKVectorDatabase &db,
                                   const Camera &camera,
                                   float tolerance) {

    assert(identifiers->size() == 4);
    StarIdentifiers pyramidIdentifiers = *identifiers; // copy with only the pyramid's high confidence stars
    Vec3 pyramidActualSpatials[4];
    for (int l = 0; l < 4; l++) {
        pyramidActualSpatials[l] = camera.CameraToSpatial(stars[pyramidIdentifiers[l].starIndex].position).Normalize();
    }

    for (int p = 0; p < (int)stars.size(); p++) {
        // ensure this star isn't in the pyramid
        bool pInPyramid = false;
        for (const StarIdentifier &id : pyramidIdentifiers) {
            if (id.starIndex == p) {
                pInPyramid = true;
                break;
            }
        }
        if (pInPyramid) {
            continue;
        }

        Vec3 pSpatial = camera.CameraToSpatial(stars[p].position).Normalize();
        float ipDist = AngleUnit(pyramidActualSpatials[0], pSpatial);
        long ipNum;
        const int16_t *ipPairs = db.FindPairsLiberal(ipDist - tolerance, ipDist + tolerance, &ipNum);
        PairDistanceInvolvingIterator pIterator(ipPairs, ipNum, pyramidIdentifiers[0].catalogIndex);

        std::vector<int16_t> pCandidates; // collect them all in the loop, at the end only identify
                                          // the star if unique
        while (pIterator.hasValue()) {
            bool ok = true;
            for (int l = 1; l < 4; l++) {
                float actualDist = AngleUnit(pSpatial, pyramidActualSpatials[l]);
                float expectedDist = AngleUnit(catalog[*pIterator].spatial,
                                               catalog[pyramidIdentifiers[l].catalogIndex].spatial);
                if (actualDist < expectedDist - tolerance || actualDist > expectedDist + tolerance) {
                    ok = false;
                }
            }
            if (ok) {
                pCandidates.push_back(*pIterator);
            }
            ++pIterator;
        }

        if (pCandidates.size() == 1) {
            identifiers->push_back(StarIdentifier(p, pCandidates[0]));
        }
        if (pCandidates.size() > 1) {
            std::cerr << "duplicate other star??" << std::endl;
            for (int16_t c : pCandidates) {
                std::cerr << catalog[c].name << std::endl;
            }
        }
    }
}

StarIdentifiers PyramidStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {

    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL || stars.size() < 4) {
        std::cerr << "Not enough stars, or database missing." << std::endl;
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    // smallest normal single-precision float is around 10^-38 so we should be all good. See
    // Analytic_Star_Pattern_Probability on the HSL wiki for details.
    float expectedMismatchesConstant = pow(numFalseStars, 4) * pow(tolerance, 5) / 2 / pow(M_PI, 2);

    // this iteration technique is described in the Pyramid paper. Briefly: i will always be the
    // lowest index, then dj and dk are how many indexes ahead the j-th star is from the i-th, and
    // k-th from the j-th. In addition, we here add some other numbers so that the pyramids are not
    // weird lines in wide FOV images. TODO: Select the starting points to ensure that the first pyramids are all within measurement tolerance.
    int numStars = (int)stars.size();
    // the idea is that the square root is about across the FOV horizontally
    int across = floor(sqrt(numStars))*2;
    int halfwayAcross = floor(sqrt(numStars)/2);
    long totalIterations = 0;

    int jMax = numStars - 3;
    for (int jIter = 0; jIter < jMax; jIter++) {
        int dj = 1+(jIter+halfwayAcross)%jMax;

        int kMax = numStars-dj-2;
        for (int kIter = 0; kIter < kMax; kIter++) {
            int dk = 1+(kIter+across)%kMax;

            int rMax = numStars-dj-dk-1;
            for (int rIter = 0; rIter < rMax; rIter++) {
                int dr = 1+(rIter+halfwayAcross)%rMax;

                int iMax = numStars-dj-dk-dr-1;
                for (int iIter = 0; iIter <= iMax; iIter++) {
                    int i = (iIter + iMax/2)%(iMax+1); // start near the center of the photo

                    // identification failure due to cutoff
                    if (++totalIterations > cutoff) {
                        std::cerr << "Cutoff reached." << std::endl;
                        return identified;
                    }

                    int j = i+dj;
                    int k = j+dk;
                    int r = k+dr;

                    assert(i!=j && j!=k && k!=r && i!=k && i!=r && j!=r);

                    // TODO: move this out of the loop?
                    Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
                    Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
                    Vec3 kSpatial = camera.CameraToSpatial(stars[k].position).Normalize();

                    float ijDist = AngleUnit(iSpatial, jSpatial);

                    float iSinInner = sin(Angle(jSpatial - iSpatial, kSpatial - iSpatial));
                    float jSinInner = sin(Angle(iSpatial - jSpatial, kSpatial - jSpatial));
                    float kSinInner = sin(Angle(iSpatial - kSpatial, jSpatial - kSpatial));

                    // if we made it this far, all 6 angles are confirmed! Now check
                    // that this match would not often occur due to chance.
                    // See Analytic_Star_Pattern_Probability on the HSL wiki for details
                    float expectedMismatches = expectedMismatchesConstant
                        * sin(ijDist)
                        / kSinInner
                        / std::max(std::max(iSinInner, jSinInner), kSinInner);

                    if (expectedMismatches > maxMismatchProbability) {
                        std::cout << "skip: mismatch prob." << std::endl;
                        continue;
                    }

                    Vec3 rSpatial = camera.CameraToSpatial(stars[r].position).Normalize();

                    // sign of determinant, to detect flipped patterns
                    bool spectralTorch = iSpatial.crossProduct(jSpatial)*kSpatial > 0;

                    float ikDist = AngleUnit(iSpatial, kSpatial);
                    float irDist = AngleUnit(iSpatial, rSpatial);
                    float jkDist = AngleUnit(jSpatial, kSpatial);
                    float jrDist = AngleUnit(jSpatial, rSpatial);
                    float krDist = AngleUnit(kSpatial, rSpatial); // TODO: we don't really need to
                                                                  // check krDist, if k has been
                                                                  // verified by i and j it's fine.

                    // we check the distances with the extra tolerance requirement to ensure that
                    // there isn't some pyramid that's just outside the database's bounds, but
                    // within measurement tolerance of the observed pyramid, since that would
                    // possibly cause a non-unique pyramid to be identified as unique.
#define _CHECK_DISTANCE(_dist) if (_dist < vectorDatabase.MinDistance() + tolerance || _dist > vectorDatabase.MaxDistance() - tolerance) { continue; }
                    _CHECK_DISTANCE(ikDist);
                    _CHECK_DISTANCE(irDist);
                    _CHECK_DISTANCE(jkDist);
                    _CHECK_DISTANCE(jrDist);
                    _CHECK_DISTANCE(krDist);
#undef _CHECK_DISTANCE

                    long ijNum, ikNum, irNum; //, jkNum, jrNum, krNum;
                    const int16_t *ijQuery = vectorDatabase.FindPairsLiberal(ijDist - tolerance, ijDist + tolerance, &ijNum);
                    const int16_t *ikQuery = vectorDatabase.FindPairsLiberal(ikDist - tolerance, ikDist + tolerance, &ikNum);
                    const int16_t *irQuery = vectorDatabase.FindPairsLiberal(irDist - tolerance, irDist + tolerance, &irNum);


                    int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
                    std::vector<bool> iSeen(catalog.size(), false);
                    for (int p = 0; p < ijNum*2; p++) {
                        int iCandidate = ijQuery[p];
                        if (iSeen[iCandidate]) {
                            continue;
                        }
                        iSeen[iCandidate] = true;

                        const Vec3 &iCandidateSpatial = catalog[iCandidate].spatial;

                        // TODO: caching these iterator results into vectors can improve
                        // performance, but at the cost of memory. It would be best to put some kind
                        // of guarantee on the memory usage, and then switch to using the iterator
                        // without caching if that memory limit is exceeded.
                        PairDistanceInvolvingIterator jIterator(ijQuery, ijNum, iCandidate);
                        PairDistanceInvolvingIterator kIterator(ikQuery, ikNum, iCandidate);
                        PairDistanceInvolvingIterator rIterator(irQuery, irNum, iCandidate);
                        std::vector<int16_t> jCandidates;
                        std::vector<int16_t> kCandidates;
                        std::vector<int16_t> rCandidates;
                        while (jIterator.hasValue()) {
                            jCandidates.push_back(*jIterator);
                            ++jIterator;
                        }
                        while (kIterator.hasValue()) {
                            kCandidates.push_back(*kIterator);
                            ++kIterator;
                        }
                        while (rIterator.hasValue()) {
                            rCandidates.push_back(*rIterator);
                            ++rIterator;
                        }
                        // TODO: break fast if any of the iterators are empty, if it's any
                        // significant performance improvement.
                        for (int16_t jCandidate : jCandidates) {
                            const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;
                            Vec3 ijCandidateCross = iCandidateSpatial.crossProduct(jCandidateSpatial);

                            for (int16_t kCandidate : kCandidates) {
                                Vec3 kCandidateSpatial = catalog[kCandidate].spatial;
                                bool candidateSpectralTorch = ijCandidateCross*kCandidateSpatial > 0;
                                // checking the spectral-ity early to fail fast
                                if (candidateSpectralTorch != spectralTorch) {
                                    continue;
                                }

                                // small optimization: We can calculate jk before iterating through r, so we will!
                                float jkCandidateDist = AngleUnit(jCandidateSpatial, kCandidateSpatial);
                                if (jkCandidateDist < jkDist - tolerance || jkCandidateDist > jkDist + tolerance) {
                                    continue;
                                }

                                // TODO: if there are no jr matches, there's no reason to
                                // continue iterating through all the other k-s. Possibly
                                // enumarete all r matches, according to ir, before this loop
                                for (int16_t rCandidate : rCandidates) {
                                    const Vec3 &rCandidateSpatial = catalog[rCandidate].spatial;
                                    float jrCandidateDist = AngleUnit(jCandidateSpatial, rCandidateSpatial);
                                    float krCandidateDist;
                                    if (jrCandidateDist < jrDist - tolerance || jrCandidateDist > jrDist + tolerance) {
                                        continue;
                                    }
                                    krCandidateDist = AngleUnit(kCandidateSpatial, rCandidateSpatial);
                                    if (krCandidateDist < krDist - tolerance || krCandidateDist > krDist + tolerance) {
                                        continue;
                                    }

                                    // we have a match!

                                    if (iMatch == -1) {
                                        iMatch = iCandidate;
                                        jMatch = jCandidate;
                                        kMatch = kCandidate;
                                        rMatch = rCandidate;
                                    } else {
                                        // uh-oh, stinky!
                                        // TODO: test duplicate detection, it's hard to cause it in the real catalog...
                                        std::cerr << "Pyramid not unique, skipping..." << std::endl;
                                        goto sensorContinue;
                                    }
                                }
                            }

                        }
                    }

                    if (iMatch != -1) {
                        printf("Matched unique pyramid!\nExpected mismatches: %e\n", expectedMismatches);
                        identified.push_back(StarIdentifier(i, iMatch));
                        identified.push_back(StarIdentifier(j, jMatch));
                        identified.push_back(StarIdentifier(k, kMatch));
                        identified.push_back(StarIdentifier(r, rMatch));

                        PyramidIdentifyRemainingStars(&identified, stars, catalog, vectorDatabase, camera, tolerance);
                        printf("Identified an additional %d stars\n", (int)identified.size() - 4);

                        return identified;
                    }

                sensorContinue:;
                }
            }
        }
    }

    std::cerr << "Tried all pyramids; none matched." << std::endl;
    return identified;
}

/* Ratio of bin size to error range, where bins are the discretization of */
/* Pattern data values to allow them to be hashed and where the error range covers */
/* a range of twice the data value's maximum error.  When a data value's error range */
/* overlaps two bins, it's replicated into both.  By linearity of expectation, */
/* the expected number of replicas of a given Pattern is: */
/* (1+1/bin_size_ratio)^num_dims, where num_dims is the number of data values */
/* stored in Patterns.  The expected ratio of Patterns with matching bins to */
/* Patterns with matching values is: (1 + bin_size_ratio)^num_dims. */
/* The bin_size_ratio represents a tradeoff between catalog size and */
/* runtime, as more replicas means a larger catalog and more mismatching results */
/* means more time spent checking for a match.  In most cases, Patterns with */
/* matching values are rare, so a larger bin_size_ratio is better.  Must be */
/* greater than 1.0 without changes to Pattern struct, recommended value of 2.0. */
#define bin_size_ratio 3.0

/* Maximum star coordinate centroiding error as fraction of maximum FOV. */
/* .001 is .1% of the max FOV or 1.414 pixels in a 1000x1000 image. */
// 1 / (1024*sqrt(2)) < .00069054
#define max_centroid_error .00069054

/* Maximum Field of View for catalog in radians.  Must exactly */
/* match the max_fov used to generate the catalog . */
#define max_fov .247

/* Maximum error in imager FOV estimate as fraction of true imager FOV. */
/* Must be greater than zero to prevent division by zero errors. */
/* max_fov*(1+max_fov_error) must be less than pi, but should be much less. */
/* Should be less than or equal to the max_fov_error used to generate the catalog. */
/* .01 for a 10 degree FOV imager covers estimates from 9.9 to 10.1 degrees. */
#define max_fov_error 0.01

/* The following values are not user defined constants and should not be changed. */
/* Maximum scaling of image caused by FOV error. */
float max_scale_factor = fmax(tan(max_fov*(1+max_fov_error)/2.0)/tan(max_fov/2.0),
                              1-tan(max_fov*(1-max_fov_error)/2.0)/tan(max_fov/2.0));

/* Largest edge error is proportional to the largest edge ratio by the image scale. */
/* For example, doubling the largest edge length doubles the scaling error. */
#define le_error_slope (max_scale_factor-1)

/* Largest edge error has a constant factor determined by the centroiding error. */
/* This value is the worst case of maximum FOV and centroid error.  It corresponds to */
/* The FOV taking on its minimum value and both centroids having maximal error inwards. */
#define le_error_offset (2*max_centroid_error/(2-max_scale_factor))


#define NUM_STARS_IN_PATTERN 4
#define max_le_length (2*sin(max_fov*(1+max_fov_error)/2.0))

StarIdentifiers TetraStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog, const Camera &camera) const {
    
    StarIdentifiers identified;

    // Tetra.c to TetraStarIdAlgo:
    // pattern_catalog --is-- catalog
    // matches         --is-- identified

    // create vector of spacial vectors of every centroided star:
    std::vector<Vec3> imageStars; // vector of spatial star vectors
    for (int i = 0; i < stars.size(); i++) {
        imageStars.push_back(camera.CameraToSpatial(stars[i].position).Normalize());
    }


    return identified;
}

/* Calculate the exponent base for logarithmic binning. */
/* Also bounds check error values. */
static double GetBase(double error_slope,
                      double error_offset){
  /* If the fixed error is non-positive, return error and exit. */
  if(error_offset <= 0){
    printf("\nNon-positive error value detected: increase error values.\n");
    exit(EXIT_FAILURE);
  }
  /* Calculate and return base of logarithmic binning function. */
  double base = (1+error_slope)/fmax(1-error_slope, 0);
  return base;
}

/* Retrieve minimum possible pre-binned value from a logarithmically binned value. */
static double LogUnbin(int bin,
                       double error_slope,
                       double error_offset){
  /* If either of the error values is infinite, 0 is the minimum value. */
  if(!isfinite(error_slope) || !isfinite(error_offset)){
    return 0;
  }
  /* Calculate base of logarithmic binning function. */
  double base = GetBase(error_slope, error_offset);
  double min_input;
  /* Slopes much smaller than the offset result in linear binning. */
  if(base <= 1+error_offset*bin_size_ratio/10.){
    min_input = bin*2*(error_slope+error_offset)*bin_size_ratio;
  }
  /* Otherwise, calculate minimum possible logarithmic pre-binned value. */
  else{
    min_input = (pow(base, bin*bin_size_ratio)-1)*error_offset/error_slope;
  }
  return min_input;
}

/* Evenly bin largest edge length values using the fact that the maximum error is linear */
/* with respect to length.  Logarithmic functions satisfy this bin size constraint. */
static int BinLargestEdge(unsigned int largest_edge,
                            int error_ratio){
  /* Convert largest edge length back to double between 0 and 1 from implicitly divided double. */
  /* Resulting value is a ratio of max_le_length, the maximum largest edge length. */
  double le_ratio = largest_edge/((1<<16)-1.0);
  /* Adjust error in le_ratio based on error_ratio between -1 and 1. */
  /* An error_ratio of -1 will give the lower bin, +1 will give the upper bin. */
  le_ratio += error_ratio*(le_ratio*le_error_slope+le_error_offset);
  /* Find and return largest edge bin using logarithmic binning function. */
  return LogUnbin(le_ratio, le_error_slope, le_error_offset);
}

/* Retrieves catalog pattern matching image Pattern.  Returns 1 if a unique */
/* match is found.  Returns 0 if multiple matches or no matches are found. */
static int get_matching_pattern(Pattern image_pattern,
                                Pattern *catalog_pattern,
                                FILE *pattern_catalog){
	/* Cache of Pattern instances from catalog. */
  static Pattern catalog_pattern_cache[pattern_cache_size];
  /* Initialize offset of the Pattern's probing sequence in the pattern cache. */
  int cache_offset = 0;
  /* Spacing between Patterns in the same probing sequence.  Grows linearly, */
  /* resulting in quadratic probing. (i.e. 0,1,3,6,10,15...) */
  int probe_step = 1;
  /* Boolean representing whether or not a catalog match has been found yet. */
  /* A match is a catalog Pattern within max_coord_error of the given image Pattern. */
  int found_match = 0;
  /* Initialize offset of the beginning of the Pattern cache in the pattern catalog. */
  uint64_t offset = hash_pattern(image_pattern);
  /* Initialize cache of catalog Patterns. */
  _fseeki64(pattern_catalog, offset*sizeof(Pattern), SEEK_SET);
  fread(catalog_pattern_cache, sizeof(Pattern), pattern_cache_size, pattern_catalog);
  /* Iterate over catalog locations in the image Pattern's probing sequence until */
  /* a catalog location without a Pattern is found, which means no matches exist. */
  /* If the probing sequence contains Patterns with the same sub-bins, iterate up */
  /* to the last one, returning success if a single matching Pattern is found, */
  /* and returning failure if no matching Patterns are found before the last one. */
  /* If two or more matching catalog Patterns are found, exit early and */
  /* return failure, as a unique identification cannot be made. */
  while(catalog_pattern_cache[cache_offset].has_pattern){
    /* Only examine catalog Patterns with the same sub-bins as the image Pattern, */
    /* as all matches will also have the same sub-bins as the image Pattern. */
    if(hash_same(image_pattern, catalog_pattern_cache[cache_offset])){
      /* Check whether the image and catalog Patterns are a match by checking if */
      /* their corresponding Features' coordinates are all within max_coord_error. */
      if(is_match(image_pattern, catalog_pattern_cache[cache_offset])){
        /* If a match has already been found previously, this must be the second */
        /* matching catalog Pattern found.  Return failure, as a unique */
        /* identification cannot be made without checking other stars. */
        /* Note that it may be worth checking all possible matches if */
        /* verifying a match is less costly than another catalog access. */
        if(found_match){
          /* Multiple matching Patterns were found.  Return failure. */
          return 0;
        }
        /* This must be the first matching catalog Pattern found.  Store it as the */
        /* output Pattern so it will be output once its uniqueness has been verified. */
        *catalog_pattern = catalog_pattern_cache[cache_offset];
        /* Set the flag indicating a matching Pattern has already been found. */
        found_match = 1;
      }
      /* If this is the last Pattern with the same sub-bins, there is no need to */
      /* keep searching until an empty catalog location is found, as any Patterns */
      /* beyond this point cannot be matches.  Exit early to save time. */
      if(catalog_pattern_cache[cache_offset].is_last){
        break;
      }
    }
    /* Advance to the next catalog location given by quadratic probing. */
    /* If cache_offset indexes beyond pattern_cache_size, return failure. */
    if(!increment_offset(pattern_catalog,
                         catalog_pattern_cache,
                         &offset,
                         &cache_offset,
                         &probe_step)){
      return 0;
    }
  }
  /* Exactly one matching Pattern was found.  Return success. */
  if(found_match){
    return 1;
  }
  /* No matching Patterns were found.  Return failure. */
  return 0;
}


bool IdentifyStars(std::vector<Vec3> &imageStars, int imageStarIds[NUM_STARS_IN_PATTERN], const Catalog &catalog, StarIdentifiers &identified) {
    int i,j;
    Pattern new_pattern;
    Pattern catalog_pattern;
    double largest_edge_length = 0.0;
    for(i = 0; i < NUM_STARS_IN_PATTERN; i++) {
        for(j = i + 1; j < NUM_STARS_IN_PATTERN; j++) {
            // dist is a "euclidian" distance evaluator for a pair of 3D vectors...
            double new_edge_length = Vec3::EucDist(imageStars[imageStarIds[i]],
                                   imageStars[imageStarIds[j]]);
            if(new_edge_length > largest_edge_length){
                largest_edge_length = new_edge_length;
                new_pattern.fixed_star_id_1 = imageStarIds[i];
                new_pattern.fixed_star_id_2 = imageStarIds[j];
            }
        }
    }
    /* Set Pattern's largest_edge_length and largest_edge_subbin.  Both are */
    /* implicitly encoded as unsigned integers with a range from 0 up to */
    /* the sine of the maximum catalog FOV to keep the catalog compact. */
    new_pattern.largest_edge = (largest_edge_length / max_le_length) * ((1 << 16) - 1);
    /* Calculate vector along x axis of Pattern's coordinate system. */
    /* The vector points from the first fixed Pattern star to the second. */
    Vec3 x_axis_vector = imageStars[new_pattern.fixed_star_id_2] - imageStars[new_pattern.fixed_star_id_1];
    /* Calculate vector along y axis of Pattern's coordinate system. */
    Vec3 y_axis_vector = imageStars[new_pattern.fixed_star_id_2].crossProduct(imageStars[new_pattern.fixed_star_id_1]);
    /* Normalize axis vectors to unit length by dividing by their magnitudes. */
    x_axis_vector.Normalize(); 
    y_axis_vector.Normalize();
    /* Use the remaining stars to initialize the Pattern's Features. */
    int feature_index = 0;
    for (i = 0; i < NUM_STARS_IN_PATTERN; i++)
    {
        /* Skip the fixed star ids, as they don't have their own Features. */
        if (imageStarIds[i] != new_pattern.fixed_star_id_1 &&
            imageStarIds[i] != new_pattern.fixed_star_id_2)
        {
            /* Set the Feature's star id to match its corresponding star. */
            new_pattern.features[feature_index].star_id = imageStarIds[i];
            /* Calculate the normalized x and y coordinates using vector projection. */
            double x = x_axis_vector * imageStars[imageStarIds[i]] / largest_edge_length;
            double y = y_axis_vector * imageStars[imageStarIds[i]] / largest_edge_length;
            /* Set Feature's coordinates by converting to implicitly divided integers. */
            new_pattern.features[feature_index].x = x * ((1 << 14) - 1);
            new_pattern.features[feature_index].y = y * ((1 << 14) - 1);
            /* Disallow 0, as rotational ambiguity correction would fail. */
            if (new_pattern.features[feature_index].x == 0)
            {
                new_pattern.features[feature_index].x = 1;
            }
            if (new_pattern.features[feature_index].y == 0)
            {
                new_pattern.features[feature_index].y = 1;
            }
            feature_index++;
        }
    }
    /* Variable encoding which 180 degree rotation will be inserted into the catalog. */
    /* A negative value means the current rotation will be inserted. */
    /* A positive value means the opposite rotation will be inserted. */
    /* A value of zero means both rotations will be inserted into the catalog. */
    int pattern_rotation;
    /* Compute largest edge bin for use in sorting Features based on x and y bins. */
    unsigned int le_bin = BinLargestEdge(new_pattern.largest_edge, 0);
    /* Helper function for sorting Features.  Sorts by x bin, then by y bin. */
    /* Returns a positive number if the first Feature has larger bin values, */
    /* returns a negative number if the second Feature has larger bin values, */
    /* and raises an error if both Features have the same bin values. */
    int CompareBins(const void *p, const void *q) {
        /* Compare the Features' x bins first, then their y bins. */
        int p_y_bin = BinY(((Feature *)p)->y, le_bin, 0);
        int q_y_bin = BinY(((Feature *)q)->y, le_bin, 0);
        int p_x_bin = BinX(((Feature *)p)->x, le_bin, p_y_bin, 0);
        int q_x_bin = BinX(((Feature *)q)->x, le_bin, q_y_bin, 0);
        /* If the x bins have different values, the y bins don't need to be computed. */
        if (p_x_bin != q_x_bin) {
            return p_x_bin - q_x_bin;
        }
        return p_y_bin - q_y_bin;
    }
    /* Sort Pattern's Features based on coordinate bins to give a unique ordering. */
    qsort(new_pattern.features, num_stars_in_pattern - 2, sizeof(Feature), compare_bins);
    /* Create a copy of the first Feature of the Pattern. */
    Feature first_feature = new_pattern.features[0];
    /* Rotate the copy by 180 degrees by taking complements of its coordinates. */
    first_feature.x = -first_feature.x;
    first_feature.y = -first_feature.y;
    /* Compare with the last Feature's bins to determine which has the largest */
    /* x bin (with y bin as a tie breaker).  This will determine the */
    /* 180 degree rotation of the Pattern's coordinate system.  The orientation */
    /* which gives the larger Feature a positive x bin value is chosen. */
    /* Put another way, the Feature furthest from the y-axis is placed on the right. */
    /* In the case that the first and last Features' bins are ambiguous after */
    /* rotating the first Feature by 180 degrees, both orientations are inserted. */
    pattern_rotation = CompareBins((void *)&first_feature,
                                    (void *)&(new_pattern.features[NUM_STARS_IN_PATTERN - 3]));
    /* If the current rotation is incorrect, rotate the Pattern by 180 degrees by taking */
    /* the complement of its Features' bin offsets and coordinates, reversing the order */
    /* of its Features, and swapping its fixed stars before inserting it into the catalog. */
    if (pattern_rotation >= 0) {
        for (i = 0; i < NUM_STARS_IN_PATTERN - 2; i++) {
            /* Take the complement of each Feature's coordinates. */
            new_pattern.features[i].x = -new_pattern.features[i].x;
            new_pattern.features[i].y = -new_pattern.features[i].y;
        }
        /* Reverse the order of the Pattern's Features by swapping across the middle. */
        for (i = 0; i < (NUM_STARS_IN_PATTERN - 2) / 2; i++) {
            Feature feature_swap = new_pattern.features[i];
            new_pattern.features[i] = new_pattern.features[NUM_STARS_IN_PATTERN - 3 - i];
            new_pattern.features[NUM_STARS_IN_PATTERN - 3 - i] = feature_swap;
        }
        /* Swap the order of the Pattern's fixed star ids and magnitudes. */
        unsigned int fixed_star_id_swap = new_pattern.fixed_star_id_1;
        new_pattern.fixed_star_id_1 = new_pattern.fixed_star_id_2;
        new_pattern.fixed_star_id_2 = fixed_star_id_swap;
    }

    /* Check cached section of catalog for Patterns matching image Pattern. */
    if (!get_matching_pattern(new_pattern, &catalog_pattern, pattern_catalog)) {
        return 0;
    }
    /* Create matching pairs of stars by corresponding fixed_star_ids and */
    /* Feature star ids between the image Pattern and catalog Pattern.  */
    // matches[0][0] = new_pattern.fixed_star_id_1;
    // matches[1][0] = new_pattern.fixed_star_id_2;
    // matches[0][1] = catalog_pattern.fixed_star_id_1;
    // matches[1][1] = catalog_pattern.fixed_star_id_2;
    for (i = 0; i < NUM_STARS_IN_PATTERN; i++) {
        // matches[i + 2][0] = new_pattern.features[i].star_id;
        // matches[i + 2][1] = catalog_pattern.features[i].star_id;
        identified.push_back(StarIdentifier(new_pattern.features[i].star_id, catalog_pattern.features[i].star_id));
    }
    return 1;
}

bool IdentifyImage(std::vector<Vec3> imageStars, const Catalog &catalog, int num_image_stars, StarIdentifiers &identified, int num_stars_selected) {
    // array of star ID's for a given pattern
    static int imageStarIds[NUM_STARS_IN_PATTERN];
    // recursively select all the image pattern stars:
    if (num_stars_selected < NUM_STARS_IN_PATTERN) {
        for (imageStarIds[num_stars_selected] = NUM_STARS_IN_PATTERN-num_stars_selected-1; 
             imageStarIds[num_stars_selected] < num_image_stars;
             imageStarIds[num_stars_selected]++) {
            
            if (IdentifyImage(imageStars, catalog, imageStarIds[num_stars_selected], identified, num_stars_selected+1)) {
                return true;
            }
        }
    } else {
        /* need to implement IdentifyStars() */
        if (IdentifyStars(imageStars, imageStarIds, catalog, identified)) {
            return true;
        }
    }
    return false;
}

}
