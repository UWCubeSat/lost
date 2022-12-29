#include "star-id.hpp"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <algorithm>
#include <fstream>  // remove after integrating database for Tetra
#include <set>
#include <utility>  // std::pair

#include "attitude-utils.hpp"
#include "databases.hpp"

namespace lost {


static bool GetCentroidCombination(std::vector<int>* const res, int pattSize, int numCentroids);


// TODO: duplicate in databases.cpp
int TetraStarIdAlgorithm::KeyToIndex(std::vector<int> key, int binFactor, int maxIndex) const {
  // key = 5-tuple of binned edge ratios
  // Outputs a row of the Pattern Catalog
  long index = 0;
  for (int i = 0; i < (int)key.size(); i++) {
    index += key[i] * std::pow(binFactor, i);
  }
  return (index * MAGIC_RAND) % maxIndex;
}

// std::vector<std::vector<int>> TetraStarIdAlgorithm::GetAtIndex(int index,
// TetraDatabase db) const{
// std::vector<std::vector<int>> TetraStarIdAlgorithm::GetAtIndex(
//     int index, std::ifstream &pattCatFile) const {

std::vector<std::vector<int>> TetraStarIdAlgorithm::GetAtIndex(int index,
                                                               const TetraDatabase &db) const {
  // Returns a list of rows (4-tuples) from the Pattern Catalog
  // Does quadratic probing

  int maxInd = catalogLength;
  // std::ifstream pattCatFile("pattCat.bin", std::ios_base::binary);

  std::vector<std::vector<int>> res;
  for (int c = 0;; c++) {
    int i = (index + c * c) % maxInd;

    // std::vector<int> tableRow = db.pattCatalog[i];
    // std::vector<int> tableRow;
    std::vector<int> tableRow = db.GetPattern(i);

    // short *row = new short[numPattStars];
    // pattCatFile.seekg(sizeof(short) * numPattStars * i, std::ios::beg);
    // pattCatFile.read((char *)row, sizeof(short) * numPattStars);

    // for (int j = 0; j < numPattStars; j++) {
    //     tableRow.push_back(int(row[j]));
    // }

    if (std::all_of(tableRow.begin(), tableRow.end(), [](int ele) { return ele == 0; })) {
    //   std::cout << c << std::endl;
      break;
    } else {
      res.push_back(tableRow);
    }
  }
  return res;
}

static bool GetCentroidCombination(std::vector<int> *const res, int pattSize, int numCentroids){
  if(numCentroids < pattSize){
    return false;
  }
  static bool firstTime = true;
  static std::vector<int> indices;
  if(firstTime){
    firstTime = false;
    indices.push_back(-1);
    for(int i = 0; i < pattSize; i++){
      indices.push_back(i);
    }
    indices.push_back(numCentroids);
    *res = std::vector<int>(indices.begin() + 1, indices.end() - 1);
    return true;
  }

  if(indices[1] >= numCentroids - pattSize){
    return false;
  }

  for(int i = 1; i <= pattSize; i++){
    indices[i]++;
    if(indices[i] < indices[i+1]){
      break;
    }
    indices[i] = indices[i-1] + 1;
  }
  *res = std::vector<int>(indices.begin() + 1, indices.end() - 1);
  return true;

}

/**
 * @brief Identifies stars using Tetra Star Identification Algorithm
 *
 * @param database Pattern catalog
 * @param stars List of star centroids detected from image
 * @param catalog Star table
 * @param camera // TODO: configure camera module for Tetra support
 * @return StarIdentifiers
 */
StarIdentifiers TetraStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars,
                                         const Catalog &catalog, const Camera &camera) const {
  // format: (centroidIndex, catalogIndex)
  StarIdentifiers result;
  std::cout << "TETRA" << std::endl;

  MultiDatabase multiDatabase(database);
  const unsigned char *databaseBuffer =
      multiDatabase.SubDatabasePointer(TetraDatabase::kMagicValue);
  if (databaseBuffer == NULL) {
    std::cerr << "Could not get to Tetra database" << std::endl;
    return result;
  }
  TetraDatabase tetraDatabase(databaseBuffer);

  int catLength = tetraDatabase.Size(); // TODO: use this, not hardcoded value in star-id.hpp


  std::vector<int> centroidIndices;
  for (int i = 0; i < (int)stars.size(); i++) {
    centroidIndices.push_back(i);
  }

  // Sort centroided stars by brightness, high to low

  std::stable_sort(centroidIndices.begin(), centroidIndices.end(),
                   [&stars](int a, int b) { return stars[a].magnitude > stars[b].magnitude; });

  // TODO: implement the generator function
  // Currently doing a naive way,  just taking the first 4 centroids
  // "Generator function" can really just be another fancy loop
  // tetra does this somewhat strangely
  // TODO: pattCheckingStars = 6, number of stars used to create possible patterns for lookup in db
  // const int pattCheckingStars = 6;



  std::vector<Star> chosenStars;
  for (int i = 0; i < numPattStars; i++) {
    chosenStars.push_back(stars[centroidIndices[i]]);
  }

  // Compute Vec3 spatial vectors for each star chosen
  // to be in the Pattern
  std::vector<Vec3> pattStarVecs;  // size==numPattStars
  for (const Star &star : chosenStars) {
    // Vec3 spatialVec = camera.CameraToSpatialFov(star.position);
    // CameraToSpatial produces different vector but also works out in the end
    // TODO: examine why the math works this way
    Vec3 spatialVec = camera.CameraToSpatial(star.position);
    pattStarVecs.push_back(spatialVec);
  }

  // Compute angle between each pair of stars chosen to be in the Pattern
  // If any angle > maxFov, then we should throw away this Pattern
  // choice, since our database will not store it
  bool angleAcceptable = true;
  for (int i = 0; i < (int)pattStarVecs.size(); i++) {
    Vec3 u = pattStarVecs[i];
    for (int j = i + 1; j < (int)pattStarVecs.size(); j++) {
      Vec3 v = pattStarVecs[j];
      float angle = Angle(u, v);  // in radians
      if (RadToDeg(angle) > maxFov) {
        angleAcceptable = false;
        break;
      }
    }
  }

  if (!angleAcceptable) {
    std::cerr << "Error: some angle is greater than maxFov" << std::endl;
    return result;
    // TODO: probably continue (try again) instead of returning, change
    // after implementing generator function
  }

  // Calculate all edge lengths in order to find value of largest edge
  // Since each Pattern consists of size==numPattStars stars, there will be
  // C(numPattStars, 2) edges For default 4-star Patterns, calculate C(4,
  // 2) = 6 edge lengths
  std::vector<float> pattEdgeLengths;  // default size = 6
  for (int i = 0; i < (int)pattStarVecs.size(); i++) {
    for (int j = i + 1; j < (int)pattStarVecs.size(); j++) {
      Vec3 diff = pattStarVecs[i] - pattStarVecs[j];
      pattEdgeLengths.push_back(diff.Magnitude());
    }
  }
  std::sort(pattEdgeLengths.begin(), pattEdgeLengths.end());
  float pattLargestEdge = pattEdgeLengths[(int)(pattEdgeLengths.size()) - 1];  // largest edge value

  // Divide each edge length by pattLargestEdge for edge ratios
  // size() = C(numPattStars, 2) - 1
  // We don't calculate ratio for largest edge, which would just be 1
  std::vector<float> pattEdgeRatios;
  for (int i = 0; i < (int)pattEdgeLengths.size() - 1; i++) {
    pattEdgeRatios.push_back(pattEdgeLengths[i] / pattLargestEdge);
  }

  // Binning step - discretize values so they can be used for hashing
  // We account for potential error in calculated values by testing hash codes
  // in a range
  std::vector<std::pair<int, int>> hcSpace;
  for (float edgeRatio : pattEdgeRatios) {
    std::pair<int, int> range;
    int lo = int((edgeRatio - pattMaxError) * numPattBins);
    lo = std::max(lo, 0);
    int hi = int((edgeRatio + pattMaxError) * numPattBins);
    hi = std::min(hi + 1, numPattBins);
    range = std::make_pair(lo, hi);
    hcSpace.push_back(range);
  }

  std::set<std::vector<int>> finalCodes;
  // TODO: if we maintain constant that numPattStars==4, we don't need to change this
  // Most star ID papers recommend setting numPattStars=4, no need to go higher (and definitely not
  // lower) Go through all the actual hash codes
  for (int a = hcSpace[0].first; a < hcSpace[0].second; a++) {
    for (int b = hcSpace[1].first; b < hcSpace[1].second; b++) {
      for (int c = hcSpace[2].first; c < hcSpace[2].second; c++) {
        for (int d = hcSpace[3].first; d < hcSpace[3].second; d++) {
          for (int e = hcSpace[4].first; e < hcSpace[4].second; e++) {
            std::vector<int> code{a, b, c, d, e};
            std::sort(code.begin(), code.end());
            finalCodes.insert(code);
          }
        }
      }
    }
  }

  for (std::vector<int> code : finalCodes) {
    int hashIndex = KeyToIndex(code, numPattBins, catalogLength);
    // Get a list of Pattern Catalog rows with hash code == hashIndex
    // One of these Patterns in the database could be a match to our constructed Pattern
    // std::vector<std::vector<int>> matches =
    //     GetAtIndex(hashIndex, pattCatFile);
    std::vector<std::vector<int>> matches = GetAtIndex(hashIndex, tetraDatabase);

    if ((int)matches.size() == 0) {
      std::cout << "Alert: matches size = 0, continuing" << std::endl;
      continue;
    }

    for (std::vector<int> matchRow : matches) {
      // std::cout << matchRow[0] << ", " << matchRow[1] << ", " << matchRow[2] << ", " <<
      // matchRow[3] << std::endl; Construct the Pattern we found in the database
      std::vector<int> catStarIDs;
      std::vector<Vec3> catStarVecs;
      for (int star : matchRow) {
        // TODO: definitely don't do this for final version, change

        // TODO: note that star IS the catalogIndex
        // Instead of scanning through entire catalog at the end (O(11 million)), just
        // scan through the 4 stars of the pattern

        // float *row = new float[starTableRowSize];
        // starTableFile.seekg(sizeof(float) * starTableRowSize * star,
        //                     std::ios::beg);
        // starTableFile.read((char *)row,
        //                    sizeof(float) * starTableRowSize);

        // Vec3 catVec(row[2], row[3], row[4]);
        // catStarIDs.push_back(row[6]);

        Vec3 catVec = catalog[star].spatial;
        catStarIDs.push_back(catalog[star].name);

        catStarVecs.push_back(catVec);
      }

      // Make the Catalog Pattern
      std::vector<float> catEdgeLengths;
      for (int i = 0; i < (int)catStarVecs.size(); i++) {
        for (int j = i + 1; j < (int)catStarVecs.size(); j++) {
          Vec3 diff = catStarVecs[i] - catStarVecs[j];
          catEdgeLengths.push_back(diff.Magnitude());
        }
      }
      std::sort(catEdgeLengths.begin(), catEdgeLengths.end());
      float catLargestEdge = catEdgeLengths[(int)(catEdgeLengths.size()) - 1];

      std::vector<float> catEdgeRatios;  // size() = C(numPattStars, 2) - 1
      for (int i = 0; i < (int)catEdgeLengths.size() - 1; i++) {
        catEdgeRatios.push_back(catEdgeLengths[i] / catLargestEdge);
      }

      bool skipMatchRow = false;
      for (int i = 0; i < (int)catEdgeRatios.size(); i++) {
        float val = catEdgeRatios[i] - pattEdgeRatios[i];
        val = std::abs(val);
        if (val > pattMaxError) {
          skipMatchRow = true;
        }
      }

      // Very common to skip here, database Pattern is NOT a match to ours
      if (skipMatchRow) {
        // std::cout << "Alert: Match Row skipped!!!" << std::endl;
        continue;
      }

      // std::cout << "fine match" << std::endl;

      Vec3 pattCentroid(0, 0, 0);
      for (Vec3 pattStarVec : pattStarVecs) {
        pattCentroid = pattCentroid + pattStarVec;
      }
      pattCentroid = pattCentroid * (1.0 / (int)pattStarVecs.size());

      std::vector<float> pattRadii;
      for (Vec3 pattStarVec : pattStarVecs) {
        pattRadii.push_back((pattStarVec - pattCentroid).Magnitude());
      }

      std::vector<int> chosenIndices =
          std::vector<int>(centroidIndices.begin(), centroidIndices.begin() + numPattStars);

      std::vector<int> sortedCentroidIndices = ArgsortVector<int>(chosenIndices, pattRadii);

      std::vector<Vec3> pattSortedVecs = ArgsortVector<Vec3>(pattStarVecs, pattRadii);

      Vec3 catCentroid(0, 0, 0);
      for (Vec3 catStarVec : catStarVecs) {
        catCentroid = catCentroid + catStarVec;
      }
      catCentroid = catCentroid * (1.0 / (int)catStarVecs.size());

      std::vector<float> catRadii;
      for (Vec3 catStarVec : catStarVecs) {
        catRadii.push_back((catStarVec - catCentroid).Magnitude());
      }

      std::vector<int> catSortedStarIDs = ArgsortVector<int>(catStarIDs, catRadii);
      std::vector<Vec3> catSortedVecs = ArgsortVector<Vec3>(catStarVecs, catRadii);

      for (int i = 0; i < numPattStars; i++) {
        int centroidIndex = sortedCentroidIndices[i];

        // segfault
        int resultStarID = catSortedStarIDs[i];

        std::cout << "Centroid Index: " << centroidIndex << ", Result StarID: " << resultStarID
                  << std::endl;

        // TODO later: take out this function, since once we actually
        // use the given catalog, we can just give the actual index
        int catalogIndex = FindCatalogStarIndex(catalog, resultStarID);

        // const CatalogStar *catStar = FindNamedStar(catalog,
        // resultStarID);

        result.push_back(StarIdentifier(centroidIndex, catalogIndex));
      }

      // NOTE: StarIdentifier wants the catalog INDEX, not the real star
      // ID

      std::cout << "SUCCESS: stars successfully matched" << std::endl;
      return result;
    }
  }

  std::cout << "FAIL" << std::endl;

  return result;
}

StarIdentifiers DummyStarIdAlgorithm::Go(const unsigned char *, const Stars &stars,
                                         const Catalog &catalog, const Camera &) const {
  StarIdentifiers result;

  unsigned int randomSeed = 123456;
  for (int i = 0; i < (int)stars.size(); i++) {
    result.push_back(StarIdentifier(i, rand_r(&randomSeed) % catalog.size()));
  }

  return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(const unsigned char *database,
                                                   const Stars &stars, const Catalog &catalog,
                                                   const Camera &camera) const {
  StarIdentifiers identified;
  MultiDatabase multiDatabase(database);
  const unsigned char *databaseBuffer =
      multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
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
        // give a greater range for min-max Query for bigger radius
        // (GreatCircleDistance)
        float lowerBoundRange = greatCircleDistance - tolerance;
        float upperBoundRange = greatCircleDistance + tolerance;
        const int16_t *upperBoundSearch;
        const int16_t *lowerBoundSearch =
            vectorDatabase.FindPairsLiberal(lowerBoundRange, upperBoundRange, &upperBoundSearch);
        // loop from lowerBoundSearch till numReturnedPairs, add one
        // vote to each star in the pairs in the datastructure
        for (const int16_t *k = lowerBoundSearch; k != upperBoundSearch; k++) {
          if ((k - lowerBoundSearch) % 2 == 0) {
            float actualAngle = AngleUnit(catalog[*k].spatial, catalog[*(k + 1)].spatial);
            assert(actualAngle <= greatCircleDistance + tolerance * 2);
            assert(actualAngle >= greatCircleDistance - tolerance * 2);
          }
          if (!votedInPair[*k] || true) {
            // if (i == 542 && *k == 9085) {
            //     printf("INC, distance %f from query %f to %f\n",
            //     greatCircleDistance,
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
    //             printf("Star %4d received %d votes.\n", catalog[i].name,
    //             votes[i]);
    //         }
    //     }
    //     printf("Debug star: Actually voted for %d with %d votes\n",
    //            catalog[indexOfMax].name, maxVotes);
    // }
    // printf("Max votes: %d\n", maxVotes);
    // starIndex = i, catalog index = indexOfMax
    StarIdentifier newStar(i, indexOfMax);
    // Set identified[i] to value of catalog index of star w most votesr
    identified.push_back(newStar);
  }
  // optimizations? N^2
  // https://www.researchgate.net/publication/3007679_Geometric_voting_algorithm_for_star_trackers
  //
  //  Do we have a metric for localization uncertainty? Star brighntess?
  // loop i from 1 through n
  std::vector<int16_t> verificationVotes(identified.size(), 0);
  for (int i = 0; i < (int)identified.size(); i++) {
    // loop j from i+1 through n
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

      // if sDist is in the range of (distance between stars in the image
      // +- R) add a vote for the match
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

/*
 * Strategies:
 * 1. For each star, enumerate all stars which have the same combination of
 * distances to some other stars, getting down to a hopefully small (<10) list
 * of candidates for each star, then do a quad-nested loop to correlate them.
 * 2. Loop through all possible stars in the catalog for star i. Then look at
 * edge ij, using this to select possible j-th stars. If ever there is not a
 * possible j-th star, continue the i-loop. When a possible ij combination is
 * found, loop through k stars according to ik. IF none are found, continue the
 * outer i loop. If some are found, check jk for each one. For each possible ijk
 * triangle,
 */

/**
 * Given a list of star pairs, finds all those pairs which involve a certain
 * star. Here "involve" means that one of the two stars in the pair is the given
 * star.
 */
class PairDistanceInvolvingIterator {
 public:
  /**
   * Create a "past-the-end" iterator.
   * If another PairDistanceInvolvingIterator is equal to this, then it is
   * done iterating.
   */
  PairDistanceInvolvingIterator() : pairs(NULL), end(NULL){};

  /**
   * The main constructor.
   * @param pairs Start of an array of star pairs
   * @param end "past-the-end" pointer for \p pairs
   * @param The catalog index that we want to be involved in the outputted
   * pairs.
   */
  PairDistanceInvolvingIterator(const int16_t *pairs, const int16_t *end, int16_t involving)
      : pairs(pairs), end(end), involving(involving) {
    assert((end - pairs) % 2 == 0);
    ForwardUntilInvolving();
  };

  // PairDistanceInvolvingIterator operator++() {
  //     PairDistanceInvolvingIterator result(*this);
  //     ++(*this);
  //     return result;
  // }

  /// Move to the next matching pair.
  PairDistanceInvolvingIterator &operator++() {
    assert(HasValue());
    pairs += 2;
    ForwardUntilInvolving();
    return *this;
  }

  /// Access the curent pair.
  int16_t operator*() const { return curValue; }

  /// Whether the iterator is currently on a value. (false if iteration is
  /// complete)
  bool HasValue() { return pairs != end; }

  // bool operator==(const PairDistanceInvolvingIterator &other) const {
  //     return ()other.pairs == pairs;
  // }

  // bool operator!=(const PairDistanceInvolvingIterator &other) const {
  //     return !(*this == other);
  // }
 private:
  const int16_t *pairs;
  const int16_t *end;
  int16_t involving;
  int16_t curValue;

  /// like postfix++, except it's a no-op if already on a valid spot.
  void ForwardUntilInvolving() {
    while (pairs != end) {
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

/// After some stars have been identified, try to idenify the rest using a
/// faster algorithm.
void PyramidIdentifyRemainingStars(StarIdentifiers *identifiers, const Stars &stars,
                                   const Catalog &catalog, const PairDistanceKVectorDatabase &db,
                                   const Camera &camera, float tolerance) {
  assert(identifiers->size() == 4);
  StarIdentifiers pyramidIdentifiers =
      *identifiers;  // copy with only the pyramid's high confidence stars
  Vec3 pyramidActualSpatials[4];
  for (int l = 0; l < 4; l++) {
    pyramidActualSpatials[l] =
        camera.CameraToSpatial(stars[pyramidIdentifiers[l].starIndex].position).Normalize();
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
    const int16_t *ipEnd;
    const int16_t *ipPairs = db.FindPairsLiberal(ipDist - tolerance, ipDist + tolerance, &ipEnd);
    PairDistanceInvolvingIterator pIterator(ipPairs, ipEnd, pyramidIdentifiers[0].catalogIndex);

    std::vector<int16_t> pCandidates;  // collect them all in the loop, at the end only
                                       // identify the star if unique
    while (pIterator.HasValue()) {
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

StarIdentifiers PyramidStarIdAlgorithm::Go(const unsigned char *database, const Stars &stars,
                                           const Catalog &catalog, const Camera &camera) const {
  StarIdentifiers identified;
  MultiDatabase multiDatabase(database);
  const unsigned char *databaseBuffer =
      multiDatabase.SubDatabasePointer(PairDistanceKVectorDatabase::kMagicValue);
  if (databaseBuffer == NULL || stars.size() < 4) {
    std::cerr << "Not enough stars, or database missing." << std::endl;
    return identified;
  }
  PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

  // smallest normal single-precision float is around 10^-38 so we should be
  // all good. See Analytic_Star_Pattern_Probability on the HSL wiki for
  // details.
  float expectedMismatchesConstant = pow(numFalseStars, 4) * pow(tolerance, 5) / 2 / pow(M_PI, 2);

  // this iteration technique is described in the Pyramid paper. Briefly: i
  // will always be the lowest index, then dj and dk are how many indexes
  // ahead the j-th star is from the i-th, and k-th from the j-th. In
  // addition, we here add some other numbers so that the pyramids are not
  // weird lines in wide FOV images. TODO: Select the starting points to
  // ensure that the first pyramids are all within measurement tolerance.
  int numStars = (int)stars.size();
  // the idea is that the square root is about across the FOV horizontally
  int across = floor(sqrt(numStars)) * 2;
  int halfwayAcross = floor(sqrt(numStars) / 2);
  long totalIterations = 0;

  int jMax = numStars - 3;
  for (int jIter = 0; jIter < jMax; jIter++) {
    int dj = 1 + (jIter + halfwayAcross) % jMax;

    int kMax = numStars - dj - 2;
    for (int kIter = 0; kIter < kMax; kIter++) {
      int dk = 1 + (kIter + across) % kMax;

      int rMax = numStars - dj - dk - 1;
      for (int rIter = 0; rIter < rMax; rIter++) {
        int dr = 1 + (rIter + halfwayAcross) % rMax;

        int iMax = numStars - dj - dk - dr - 1;
        for (int iIter = 0; iIter <= iMax; iIter++) {
          int i = (iIter + iMax / 2) % (iMax + 1);  // start near the center of the photo

          // identification failure due to cutoff
          if (++totalIterations > cutoff) {
            std::cerr << "Cutoff reached." << std::endl;
            return identified;
          }

          int j = i + dj;
          int k = j + dk;
          int r = k + dr;

          assert(i != j && j != k && k != r && i != k && i != r && j != r);

          // TODO: move this out of the loop?
          Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
          Vec3 jSpatial = camera.CameraToSpatial(stars[j].position).Normalize();
          Vec3 kSpatial = camera.CameraToSpatial(stars[k].position).Normalize();

          float ijDist = AngleUnit(iSpatial, jSpatial);

          float iSinInner = sin(Angle(jSpatial - iSpatial, kSpatial - iSpatial));
          float jSinInner = sin(Angle(iSpatial - jSpatial, kSpatial - jSpatial));
          float kSinInner = sin(Angle(iSpatial - kSpatial, jSpatial - kSpatial));

          // if we made it this far, all 6 angles are confirmed! Now
          // check that this match would not often occur due to
          // chance. See Analytic_Star_Pattern_Probability on the HSL
          // wiki for details
          float expectedMismatches = expectedMismatchesConstant * sin(ijDist) / kSinInner /
                                     std::max(std::max(iSinInner, jSinInner), kSinInner);

          if (expectedMismatches > maxMismatchProbability) {
            std::cout << "skip: mismatch prob." << std::endl;
            continue;
          }

          Vec3 rSpatial = camera.CameraToSpatial(stars[r].position).Normalize();

          // sign of determinant, to detect flipped patterns
          bool spectralTorch = iSpatial.CrossProduct(jSpatial) * kSpatial > 0;

          float ikDist = AngleUnit(iSpatial, kSpatial);
          float irDist = AngleUnit(iSpatial, rSpatial);
          float jkDist = AngleUnit(jSpatial, kSpatial);
          float jrDist = AngleUnit(jSpatial, rSpatial);
          float krDist = AngleUnit(kSpatial, rSpatial);  // TODO: we don't really need to
                                                         // check krDist, if k has been
                                                         // verified by i and j it's fine.

          // we check the distances with the extra tolerance
          // requirement to ensure that there isn't some pyramid
          // that's just outside the database's bounds, but within
          // measurement tolerance of the observed pyramid, since that
          // would possibly cause a non-unique pyramid to be
          // identified as unique.
#define _CHECK_DISTANCE(_dist)                            \
  if (_dist < vectorDatabase.MinDistance() + tolerance || \
      _dist > vectorDatabase.MaxDistance() - tolerance) { \
    continue;                                             \
  }
          _CHECK_DISTANCE(ikDist);
          _CHECK_DISTANCE(irDist);
          _CHECK_DISTANCE(jkDist);
          _CHECK_DISTANCE(jrDist);
          _CHECK_DISTANCE(krDist);
#undef _CHECK_DISTANCE

          const int16_t *ijEnd, *ikEnd, *irEnd;
          const int16_t *const ijQuery =
              vectorDatabase.FindPairsLiberal(ijDist - tolerance, ijDist + tolerance, &ijEnd);
          const int16_t *const ikQuery =
              vectorDatabase.FindPairsLiberal(ikDist - tolerance, ikDist + tolerance, &ikEnd);
          const int16_t *const irQuery =
              vectorDatabase.FindPairsLiberal(irDist - tolerance, irDist + tolerance, &irEnd);

          int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
          std::vector<bool> iSeen(catalog.size(), false);
          for (const int16_t *iCandidateQuery = ijQuery; iCandidateQuery != ijEnd;
               iCandidateQuery++) {
            int iCandidate = *iCandidateQuery;
            if (iSeen[iCandidate]) {
              continue;
            }
            iSeen[iCandidate] = true;

            const Vec3 &iCandidateSpatial = catalog[iCandidate].spatial;

            // TODO: caching these iterator results into vectors can
            // improve performance, but at the cost of memory. It
            // would be best to put some kind of guarantee on the
            // memory usage, and then switch to using the iterator
            // without caching if that memory limit is exceeded.
            PairDistanceInvolvingIterator jIterator(ijQuery, ijEnd, iCandidate);
            PairDistanceInvolvingIterator kIterator(ikQuery, ikEnd, iCandidate);
            PairDistanceInvolvingIterator rIterator(irQuery, irEnd, iCandidate);
            std::vector<int16_t> jCandidates;
            std::vector<int16_t> kCandidates;
            std::vector<int16_t> rCandidates;
            while (jIterator.HasValue()) {
              jCandidates.push_back(*jIterator);
              ++jIterator;
            }
            while (kIterator.HasValue()) {
              kCandidates.push_back(*kIterator);
              ++kIterator;
            }
            while (rIterator.HasValue()) {
              rCandidates.push_back(*rIterator);
              ++rIterator;
            }
            // TODO: break fast if any of the iterators are empty,
            // if it's any significant performance improvement.
            for (int16_t jCandidate : jCandidates) {
              const Vec3 &jCandidateSpatial = catalog[jCandidate].spatial;
              Vec3 ijCandidateCross = iCandidateSpatial.CrossProduct(jCandidateSpatial);

              for (int16_t kCandidate : kCandidates) {
                Vec3 kCandidateSpatial = catalog[kCandidate].spatial;
                bool candidateSpectralTorch = ijCandidateCross * kCandidateSpatial > 0;
                // checking the spectral-ity early to fail fast
                if (candidateSpectralTorch != spectralTorch) {
                  continue;
                }

                // small optimization: We can calculate jk
                // before iterating through r, so we will!
                float jkCandidateDist = AngleUnit(jCandidateSpatial, kCandidateSpatial);
                if (jkCandidateDist < jkDist - tolerance || jkCandidateDist > jkDist + tolerance) {
                  continue;
                }

                // TODO: if there are no jr matches, there's no
                // reason to continue iterating through all the
                // other k-s. Possibly enumarete all r matches,
                // according to ir, before this loop
                for (int16_t rCandidate : rCandidates) {
                  const Vec3 &rCandidateSpatial = catalog[rCandidate].spatial;
                  float jrCandidateDist = AngleUnit(jCandidateSpatial, rCandidateSpatial);
                  float krCandidateDist;
                  if (jrCandidateDist < jrDist - tolerance ||
                      jrCandidateDist > jrDist + tolerance) {
                    continue;
                  }
                  krCandidateDist = AngleUnit(kCandidateSpatial, rCandidateSpatial);
                  if (krCandidateDist < krDist - tolerance ||
                      krCandidateDist > krDist + tolerance) {
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
                    // TODO: test duplicate detection, it's
                    // hard to cause it in the real
                    // catalog...
                    std::cerr << "Pyramid not unique, skipping..." << std::endl;
                    goto sensorContinue;
                  }
                }
              }
            }
          }

          if (iMatch != -1) {
            printf(
                "Matched unique pyramid!\nExpected mismatches: "
                "%e\n",
                expectedMismatches);
            identified.push_back(StarIdentifier(i, iMatch));
            identified.push_back(StarIdentifier(j, jMatch));
            identified.push_back(StarIdentifier(k, kMatch));
            identified.push_back(StarIdentifier(r, rMatch));

            PyramidIdentifyRemainingStars(&identified, stars, catalog, vectorDatabase, camera,
                                          tolerance);
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

}  // namespace lost
