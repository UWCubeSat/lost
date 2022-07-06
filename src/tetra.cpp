#include <iostream>
#include <math.h>

using namespace std;

// Number of stars in stars file
const int numStars = 5904;
// Number of pattern locations cached when catalog accessed
const int pattCacheSize = 16;
// Max distance matching pattern can be from original offset in catalog
const int maxProbeDepth = 4278;
// Number of stars in a Pattern
const int numPattStars = 4;
// Number of Patterns in Tetra catalog
const int numCatalogPatts = 770708495;
// Max FOV for catalog in radians
const double maxFov = 0.247;
// Ratio of bin size to error range
// When a data value's error range overlaps 2 bins, it is replicated into both bins
// By LOE, expected number of replicas of a given Pattern is
// (1+1/bin_size_ratio)^k, where k = number of data values in a Pattern
// Tradeoff: lower ratio means more replicas => larger catalog
// higher ratio means more collisions/mismatches => more time needed to check for a match
const double binSizeRatio = 3.0;
// Max centroiding error as fraction of maxFov TODO
const double maxCentroidError = .00069054;
// Max error in imager FOV estimate as fraction of true imager FOV TODO: what?
const double maxFovError = 0.01;
// Max number of stars to process TODO: process?
const int maxStars = 12;
// Max number of stars per image
const int maxStarsPerImg = 25;
// Number of test images
const int numImages = 1000;

// Max scaling of image caused by FOV error
const double maxScaleFactor = fmax(tan(maxFov * (1 + maxFovError) / 2.0) / tan(maxFov / 2.0),
                                   1 - tan(maxFov * (1 - maxFovError) / 2.0) / tan(maxFov / 2.0));
// Largest edge error
// TODO: don't know how these are calculated or what they are
const double leErrorSlope = maxScaleFactor - 1;
const double leErrorOffset = 2*maxCentroidError / (2-maxScaleFactor);
// Largest possible largest edge length given maxFov
const double maxLELength = 2 * sin(maxFov*(1+maxFovError) / 2.0);

class Vec3{
public:
    Vec3();
    Vec3(double x, double y, double z): x(x), y(y), z(z) {}

    Vec3 operator- (const Vec3 &other) const;
    double dot (const Vec3 &other) const;

    double mag() const;
    Vec3 normalize() const;
    double dist(const Vec3 &other) const;
    Vec3 cross(const Vec3 &other) const;

    double x;
    double y;
    double z;
};

Vec3::Vec3(){}

/**
 * Returns the 3D vector: this - other
 * This member function is const- will not modify this state
 * @param other
 * @return 3D vector = this - other
 */
Vec3 Vec3::operator- (const Vec3 &other) const{
    return {x-other.x, y-other.y, z-other.y};
}

/**
 * Returns magnitude of this
 * @return magnitude of this 3D vector
 */
double Vec3::mag() const{
    return sqrt(x*x + y*y + z*z);
}

/**
 * Returns unit vector in direction of this
 * @return 3D unit vector in direction of this
 */
Vec3 Vec3::normalize() const{
    double m = mag();
    return {x/m, y/m, z/m};
}

/**
 * Returns Euclidean distance between this and other
 * i.e. returns magnitude of this-other
 * @param other
 * @return
 */
double Vec3::dist (const Vec3 &other) const{
    Vec3 res = *this - other;
    return res.mag();
}

/**
 * Returns dot product of this and other
 * @param other
 * @return
 */
double Vec3::dot (const Vec3 &other) const{
    return x*other.x + y*other.y + z*other.z;
}

/**
 * Returns cross product of this and other
 * @param other
 * @return
 */
Vec3 Vec3::cross(const Vec3 &other) const{
    Vec3 res;
    res.x = y * other.z - other.y * z;
    res.y = other.x * z - x * other.z;
    res.z = x * other.y - other.x * y;
    return res;
}

//////////////////// STRUCTS //////////////////////
////                                        ///////
///////////////////////////////////////////////////

/**
 * Feature represents a star that is NOT part of the largest edge in a Pattern
 */
struct Feature{
    // Normalized x coordinate
    int x;
    // Normalized y coordinate
    int y;
    unsigned int xBinOffset; // TODO: what are these? Why do we need both x and y bin offset
    unsigned int yBinOffset;
    // ID of the Feature's star
    unsigned int starID; // TODO: set to weird thing in identifyStars()
};

/**
 * Pattern represents the star pattern and associated coordinate system
 * This is what is stored in a given catalog position
 * TODO: something about 180 degree rotational ambiguity
 */
struct Pattern{
    // Features represent the stars NOT part of the largest edge
    // stored in order of increasing x bin, then y bin
    Feature features[numPattStars - 2];
    // Length of largest edge in Pattern
    unsigned int leLength;
    unsigned int leBinOffset;
    // ID of first star that forms the largest edge
    unsigned int leStarID1;
    // ID of second star that forms the largest edge
    unsigned int leStarID2;
    // Indicate whether this catalog position contains the last matching
    // pattern in the catalog
    // If so, no need to probe further
    unsigned int isLast;
};

/**
 * Check if a catalog position contains a valid Pattern
 * @param p
 * @return
 */
bool hasPattern(Pattern p){
    return p.leLength > 0;
}

/**
 * Represents a star stored in the stars array
 * NOTE: this is different from a centroid star
 */
struct Star{
    // 3D unit vector pointing to star from center of celestial sphere TODO
    Vec3 vec;
    // Magnitude of the star TODO: in what units?
    double mag;
    // Star ID, i.e. Hipparchos number
    unsigned int starID;
};

/////////////////////////////////////////////////////////////////////////

/**
 * Compute exponent base for logarithmic binning
 * TODO: what is logarithmic binning? What do the params mean
 * @param errorSlope
 * @param errorOffset
 * @return
 */
double getLogBinningBase (double errorSlope, double errorOffset){
    if(errorOffset <= 0){
        cout << "Error: non-positive error value detected" << endl;
        exit(EXIT_FAILURE);
    }
    // Calculate base of logarithmic binning function
    // fmax(a, b) returns the larger of 2 floating point arguments
    double base = (1 + errorSlope) / fmax(1-errorSlope, 0);
    return base;
}

/**
 * Calculates log bin TODO: how does this work, whhat does it do
 * @param input
 * @param errorSlope
 * @param errorOffset
 * @return
 */
int logBin(double input, double errorSlope, double errorOffset){
    // If either errorSlope or errorOffset is infinite, return 0 bin
    if(isinf(errorSlope) || isinf(errorOffset)){
        return 0;
    }
    // Calculate base of logarithmic binning function
    double base = getLogBinningBase(errorSlope, errorOffset);
    int bin;
    // We do linear binning if slope << offset
    if(base <= 1 + errorOffset * binSizeRatio / 10.0){
        bin = input / (2* (errorSlope + errorOffset) * binSizeRatio);
    }else{
        bin = (log(input * errorSlope / errorOffset + 1) / log(base)) / binSizeRatio;
    }
    return bin;
}

/**
 * Compute min possible input value given its log bin
 * @param bin
 * @param errorSlope
 * @param errorOffset
 * @return
 */
double logUnbin(int bin, double errorSlope, double errorOffset){
    if(isinf(errorSlope) || isinf(errorOffset)){
        return 0;
    }
    double base = getLogBinningBase(errorSlope, errorOffset);
    double minInput;
    if(base <= 1 + errorOffset * binSizeRatio / 10.0){
        minInput = bin * 2 * (errorSlope + errorOffset) * binSizeRatio;
    }else{
        // pow(base, pwr) returns base^pwr
        minInput = (pow(base, bin*binSizeRatio) - 1) * errorOffset / errorSlope;
    }
    return minInput;
}

/**
 * Bin largest edge length- TODO: math?
 * @param largestEdgeLen
 * @param errorRatio
 * @return
 */
int binLargestEdge(unsigned int largestEdgeLen, int errorRatio){
    //TODO: what is this math?
    double leRatio = largestEdgeLen / ((1 << 16) - 1.0);
    leRatio += errorRatio * (leRatio * leErrorSlope + leErrorOffset);
    return logBin(leRatio, leErrorSlope, leErrorOffset);
}

/**
 * Return min largest edge ratio within the bin - TODO: math?
 * @param bin
 * @return
 */
double unbinLargestEdge(unsigned int bin){
    double minLERatio = logUnbin(bin, leErrorSlope, leErrorOffset);
    return minLERatio;
}

/**
 * Return y coordinate bin
 * @param y
 * @param leBin
 * @param errorRatio
 * @return
 */
int binYCoord(int y, unsigned int leBin, int errorRatio){
    double minLERatio = unbinLargestEdge(leBin);
    double errorConst = leErrorOffset / (2-maxScaleFactor);
    double errorSlope = errorConst / fmax(minLERatio-errorConst, 0);
    double errorOffset = errorSlope;

    double yRatio = y / ((1<<14) - 1.0);
    // copysign(a, b) returns a floating point value with magnitude of a & sign of b
    // fabs(x) returns absolute value of a floating point x
    yRatio += errorRatio * copysign(fabs(yRatio) * errorSlope + errorOffset, yRatio);
    int bin = logBin(fabs(yRatio), errorSlope, errorOffset);
    if(yRatio < 0){
        // ~ is bitwise NOT, inverts all bits in binary number
        bin = ~bin;
    }
    return bin;
}

/**
 * Return y coord bin into max y coordinate ratio (absolute value) within the bin
 * TODO
 * @param bin
 * @param leBin
 * @return
 */
double unbinYCoord(int bin, unsigned int leBin){
    double minLERatio = unbinLargestEdge(leBin);
    double errorConst = leErrorOffset / (2-maxScaleFactor);
    double errorSlope = errorConst / fmax(minLERatio - errorConst, 0);
    double errorOffset = errorSlope;

    double maxYRatio = logUnbin(bin >= 0 ? bin+1 : (~bin)+1, errorSlope, errorOffset);
    return maxYRatio;
}

/**
 * Bin x coordinate using y coordinate and largest edge bins
 * @param x
 * @param leBin
 * @param yBin
 * @param errorRatio
 * @return
 */
int binXCoord(int x, unsigned int leBin, int yBin, int errorRatio){
    double minLERatio = unbinLargestEdge(leBin);
    double maxYRatio = unbinYCoord(yBin, leBin);
    double errorConst = leErrorOffset / (2-maxScaleFactor);

    double errorSlope = errorConst / fmax(minLERatio-errorConst, 0);
    double errorOffset = errorSlope * (1+ 2*sqrt((1.0/4) + maxYRatio * maxYRatio)) / 2;
    double xRatio = x / ((1 << 14) - 1.0);
    xRatio += errorRatio * copysign(fabs(xRatio) * errorSlope + errorOffset, xRatio);

    int bin = logBin(fabs(xRatio), errorSlope, errorOffset);
    if(xRatio < 0){
        bin = ~bin;
    }
    return bin;
}

/**
 * Hash function TODO
 * @param oldHash
 * @param key
 * @return
 */
uint64_t hashInt(uint64_t oldHash, uint64_t key){
    key = key * 11400714819323198549ULL;
    return oldHash ^ (oldHash >> 13) ^ key;
}

/**
 * Hash function, takes a Pattern and produces a catalog position/index
 * based on the Pattern's bins
 * NOTE: this function just computes the hash value, doesn't put anything in
 * the actual catalog
 * @param patt
 * @return
 */
uint64_t hashPattern(Pattern patt){
    // init hash value to largest edge bin
    unsigned int leBin = binLargestEdge(patt.leLength, 0);
    uint64_t hash = hashInt(0, leBin);

    // Update hash using each Feature's x and y bins
    for(int i = 0; i < numPattStars-2; i++){
        int yBin = binYCoord(patt.features[i].y, leBin, 0);
        hash = hashInt(hash, yBin + (1 << 31));
        int xBin = binXCoord(patt.features[i].x, leBin, yBin, 0);
        hash = hashInt(hash, xBin + (1 << 31));
    }
    // Could result in collision
    return hash % numCatalogPatts;
}

/**
 * Checks if newPattern and catalogPattern have same bin pairs (LE, x, y)
 * in case of a collision
 * Kind of a pre-check for isMatch()
 * @param newPattern Pattern created from image
 * @param catalogPattern Pattern stored in catalog
 * @return 1 if all bin pairs are the same; 0 otherwise
 */
int checkSameBins(Pattern newPattern, Pattern catPattern){
    unsigned int newLEBin = binLargestEdge(newPattern.leLength, 0);
    unsigned int catLEBin = binLargestEdge(catPattern.leLength, 2*catPattern.leBinOffset - 1);
    if(newLEBin != catLEBin){
        return 0; // Largest edge bins don't match
    }

    for(int i = 0; i < numPattStars-2; i++){
        Feature newFeature = newPattern.features[i];
        Feature catFeature = catPattern.features[i];
        int newYBin = binYCoord(newFeature.y, newLEBin, 0);
        int catYBin = binYCoord(catFeature.y, catLEBin, 2*catFeature.yBinOffset - 1);

        if(newYBin != catYBin){
            return 0; // y bins don't match => collision
        }

        int newXBin = binXCoord(newFeature.x, newLEBin, newYBin, 0);
        int catXBin = binXCoord(catFeature.x, catLEBin, catYBin, 2*catFeature.xBinOffset - 1);

        if(newXBin != catXBin){
            return 0; // x bins don't match => collision
        }
    }

    return 1;
}

/**
 * Check whether two Patterns match by checking if their x, y coordinates and
 * largest edge lengths "match"- TODO: match check is interesting
 * @param newPattern Pattern created from image
 * @param catPattern Pattern stored in catalog
 * @return 1 if the Patterns match; 0 otherwise
 */
int isMatch(Pattern newPattern, Pattern catPattern){

    // Check that image Pattern's LE length is "within range"
    // of the catalog Pattern
    double newLERatio = newPattern.leLength / ((1 << 16) - 1.0);
    double catLERatio = catPattern.leLength / ((1 << 16) - 1.0);
    double maxLEError = catLERatio * leErrorSlope + leErrorOffset;
    if(fabs(newLERatio - catLERatio) > maxLEError){
        return 0;
    }

    double coordErrorConst = leErrorOffset / (2-maxScaleFactor);
    double coordErrorSlope = coordErrorConst / fmax(newLERatio-coordErrorConst, 0);
    double coordErrorOffsetY = coordErrorSlope;

    for(int i = 0; i < numPattStars-2; i++){
        double newYCoord = newPattern.features[i].y / ((1 << 14) - 1.0);
        double catYCoord = catPattern.features[i].y / ((1 << 14) - 1.0);
        double maxYError = fabs(catYCoord) * coordErrorSlope + coordErrorOffsetY;
        if(fabs(newYCoord - catYCoord) > maxYError){
            return 0;
        }
    }

    unsigned int catLEBin = binLargestEdge(catPattern.leLength, 2*catPattern.leBinOffset - 1);
    for(int i = 0; i < numPattStars-2; i++){
        int catYBin = binYCoord(catPattern.features[i].y, catLEBin,
                                2*catPattern.features[i].yBinOffset - 1);
        double maxYRatio = unbinYCoord(catYBin, catLEBin);
        double coordErrorOffsetX = coordErrorSlope * (1 + 2*sqrt((1.0/4) + maxYRatio * maxYRatio)) / 2;

        double newXCoord = newPattern.features[i].x / ((1 << 14) - 1.0);
        double catXCoord = catPattern.features[i].x / ((1 << 14) - 1.0);
        double maxXError = fabs(catXCoord) * coordErrorSlope + coordErrorOffsetX;
        if(fabs(newXCoord - catXCoord) > maxXError){
            return 0;
        }
    }

    return 1;
}

/**
 * Update where we are in the Pattern cache, possibly update the cache itself
 * @param pattCatalog Catalog storing all Patterns
 * @param catCache Cache that stores catalog Patterns, used in get_matching_pattern()
 * @param offset Where the cache starts in the CATALOG
 * @param cacheOffset Where we are in the CACHE
 * @param probeStep How far we step in the cache to find next possible match
 * @return 0 if probe goes out of probe founds; else return 1
 */
int incrementOffset(FILE *pattCatalog, Pattern catCache[pattCacheSize],
                    uint64_t *offset, int* cacheOffset, int* probeStep){

    if(((*probeStep) * (*probeStep + 1)) / 2 > maxProbeDepth){
        //TODO: why this calculation, n(n+1)/2
        return 0;
    }
    *cacheOffset += *probeStep;
    // If we probe outside our current cache, update cache to next probe offset
    if(*cacheOffset >= pattCacheSize){
        // Update offset (where we are in the catalog)
        *offset += *cacheOffset;
        // Reset our offset in the cache itself
        *cacheOffset = 0;
        // Moves file pointer {offset} bytes from {origin}
        // SEEK_SET = 0 = beginning of file
        // TODO: change back to fseeko64?
        fseeko64(pattCatalog, *offset * sizeof(Pattern), SEEK_SET);
        // TODO: right now, sizeof(Pattern) = 60
        // Reads an array of {pattCacheSize=16} Patterns from catalog
        // and stores them in our cache
        fread(catCache, sizeof(Pattern), pattCacheSize, pattCatalog);
    }
    *probeStep += 1;
    return 1; // probe stayed within probe bounds
}

/**
 * Looks through catalog to see if there is a catalog Pattern
 * matching the Pattern we constructed from our image
 *
 * Arguably the most important function in Tetra
 * We use a Pattern cache that stores a section of the catalog to look for matches
 *
 * @param imgPattern Our constructed Pattern that we want to find a match for
 * @param catalogPattern Output = pattern in catalog that matches ours
 * @param pattCatalog Catalog storing all Patterns
 * @return 1 if a UNIQUE match is found; else return 0 if multiple or no matches found
 */
int getMatchingPattern(Pattern imgPattern, Pattern *catalogPattern, FILE *pattCatalog){
    // TODO: why make this static
    // Cache of Patterns from the catalog
    static Pattern catCache[pattCacheSize];
    // Explore cache from beginning
    int cacheOffset = 0;
    // This grows linearly (0, 1, 3, 6, 10, 15, ...) which results in quadratic probing
    int probeStep = 1;
    // Track whether a matching catalog Pattern has been found yet
    bool foundMatch = false;
    // Initialize beginning of cache in catalog to hash of our imgPattern
    uint64_t offset = hashPattern(imgPattern);
    // Start our catalog pointer at offset * sizeof(Pattern)
    // TODO: check fseek
    fseeko64(pattCatalog, offset * sizeof(Pattern), SEEK_SET);
    // Fill our cache
    fread(catCache, sizeof(Pattern), pattCacheSize, pattCatalog);

    // // Perform quadratic probing through the catalog
    // while(hasPattern(catCache[cacheOffset]){
    //     Pattern catPattern = catCache[cacheOffset];




    // }
    return 2;
}




int main(){
    Pattern p;
    Pattern q;
    cout << sizeof(Pattern) << endl; // 60 bytes

    return 0;
}