#include "tetra-db.hpp"

int Tetra::KeyToIndex(std::vector<int> key, int binFactor, int maxIndex) {
    const long long MAGIC_RAND = 2654435761;
    long index = 0;
    for (int i = 0; i < (int)key.size(); i++) {
        index += key[i] * std::pow(binFactor, i);
    }
    return (index * MAGIC_RAND) % maxIndex;
}

void Tetra::GenerateDatabase(float maxFov, std::string &starTableFName,
                             std::string &catalogFName, std::string &propsFName,
                             int pattStarsPerFOV, int verificationStarsPerFOV,
                             float starMaxMagnitude, float starMinSeparation,
                             float pattMaxError) {
    maxFov = maxFov * (float)M_PI / 180;
    const short pattSize = 4;
    const short pattBins = 25;
    const time_t now = time(nullptr);
    const tm *utc = gmtime(&now);
    const int currYear = utc->tm_year + 1900;

    const int headerLength = 28;
    int numEntries = 9110;
    std::cout << "Loading bsc5 catalog with " << numEntries << " star entries"
              << std::endl;

    std::vector<StarTableEntry> starTable(numEntries);

    std::ifstream bsc5File("BSC5", std::ios::in | std::ios::binary);
    if (!bsc5File.is_open()) {
        std::cout << "Error: Could not open bsc5 file" << std::endl;
        return;
    }
    bsc5File.seekg(headerLength, std::ios::beg);

    std::vector<Bsc5DataEntry> reader(numEntries);
    for (int i = 0; i < numEntries; i++) {
        bsc5File.read((char *)&reader[i].id, sizeof(reader[i].id));
        bsc5File.read((char *)&reader[i].ra_1950, sizeof(reader[i].ra_1950));
        bsc5File.read((char *)&reader[i].dec_1950, sizeof(reader[i].dec_1950));
        bsc5File.read((char *)&reader[i].type, sizeof(reader[i].type));
        bsc5File.read((char *)&reader[i].magnitude,
                      sizeof(reader[i].magnitude));
        bsc5File.read((char *)&reader[i].ra_pm, sizeof(reader[i].ra_pm));
        bsc5File.read((char *)&reader[i].dec_pm, sizeof(reader[i].dec_pm));

        if (!bsc5File.good()) {
            std::cout << "Error: Could not read bsc5 file" << std::endl;
            std::cout << "Ended at entry " << i << std::endl;
            bsc5File.close();
            return;
        }
    }

    bsc5File.close();
    std::cout << "Successfully read bsc5 file" << std::endl;

    // Read into star table
    for (int i = 0; i < numEntries; i++) {
        const float mag = (float)reader[i].magnitude / 100;
        const float id = (float)reader[i].id;
        if (mag <= starMaxMagnitude) {
            float ra = (float)reader[i].ra_1950 +
                       reader[i].ra_pm * ((float)currYear - 1950);
            float dec = (float)reader[i].dec_1950 +
                        reader[i].dec_pm * ((float)currYear - 1950);

            starTable[i].id = id;
            starTable[i].ra = ra;
            starTable[i].dec = dec;
            starTable[i].magnitude = mag;
        }
    }

    // correct
    std::vector<bool> kept(numEntries);
    int keptCount = 0;
    for (int i = 0; i < numEntries; i++) {
        if (starTable[i].ra != 0 || starTable[i].dec != 0) {
            kept[i] = true;
            keptCount++;
        }
    }
    std::vector<StarTableEntry> newStarTable;
    for (int i = 0; i < numEntries; i++) {
        if (kept[i]) {
            newStarTable.push_back(starTable[i]);
        }
    }
    starTable = newStarTable;
    numEntries = keptCount;
    std::cout << "Loaded " << numEntries << " stars with magnitude below "
              << starMaxMagnitude << std::endl
              << std::endl;

    // TODO: sort by magnitude (low to high), not stable
    // This implies brighter stars come first
    // stable_sort produces same result as py
    // TODO: question now is, does stable work? - YES

    // TODO: presort in some other place
    // probably not O(1) after all this stuff
    std::stable_sort(starTable.begin(), starTable.end(),
                     [](const StarTableEntry &a, const StarTableEntry &b) {
                         return a.magnitude < b.magnitude;
                     });


    // Calculate star direction vectors
    for (int i = 0; i < numEntries; i++) {
        starTable[i].x = std::cos(starTable[i].ra) * std::cos(starTable[i].dec);
        starTable[i].y = std::sin(starTable[i].ra) * std::cos(starTable[i].dec);
        starTable[i].z = std::sin(starTable[i].dec);
    }

    // NOTE: BscParse() and NarrowCatalog() combined take care of all steps prior to this point

    // NOTE: beginning of new code


    int keepForPattCount = 1;
    std::vector<bool> keepForPatterns(numEntries);
    std::vector<bool> keepForVerifying(numEntries);
    keepForPatterns[0] = true;
    keepForVerifying[0] = true;
    // allStarVectors is 3 x 9050
    for (int ind = 1; ind < numEntries; ind++) {
        std::vector<float> vec{starTable[ind].x, starTable[ind].y,
                               starTable[ind].z};

        std::vector<float> angsPatterns;
        for (int j = 0; j < numEntries; j++) {
            if (keepForPatterns[j]) {
                float dotProd = vec[0] * starTable[j].x +
                                vec[1] * starTable[j].y +
                                vec[2] * starTable[j].z;
                angsPatterns.push_back(dotProd);
            }
        }

        std::vector<float> angsVerifying;
        for (int j = 0; j < numEntries; j++) {
            if (keepForVerifying[j]) {
                float dotProd = vec[0] * starTable[j].x +
                                vec[1] * starTable[j].y +
                                vec[2] * starTable[j].z;
                angsVerifying.push_back(dotProd);
            }
        }

        bool angsPatternsOK = true;
        for (float angPatt : angsPatterns) {
            if (angPatt >= std::cos(starMinSeparation * (float)M_PI / 180)) {
                angsPatternsOK = false;
                break;
            }
        }
        if (std::isnan(starMinSeparation) || angsPatternsOK) {
            int numStarsInFOV = 0;
            for (float angPatt : angsPatterns) {
                if (angPatt > std::cos(maxFov / 2)) numStarsInFOV++;
            }
            if (numStarsInFOV < pattStarsPerFOV) {
                keepForPatterns[ind] = true;
                keepForVerifying[ind] = true;
                keepForPattCount++;
            }
        }

        bool angsVerifyingOK = true;
        for (float angVer : angsVerifying) {
            if (angVer >= std::cos(starMinSeparation * (float)M_PI / 180)) {
                angsVerifyingOK = false;
                break;
            }
        }
        if (std::isnan(starMinSeparation) || angsVerifyingOK) {
            int numStarsInFOV = 0;
            for (float angVer : angsVerifying) {
                if (angVer > std::cos(maxFov / 2)) numStarsInFOV++;
            }
            if (numStarsInFOV < verificationStarsPerFOV) {
                keepForVerifying[ind] = true;
            }
        }
    }

    std::vector<StarTableEntry> finalStarTable;
    for (int i = 0; i < (int)keepForVerifying.size(); i++) {
        if (keepForVerifying[i]) {
            finalStarTable.push_back(starTable[i]);
        }
    }
    starTable = finalStarTable;

    // Good!
    // Star Table has 7099 entries
    std::vector<short> pattStars;
    short cumulativeSum = -1;
    for (int i = 0; i < (int)keepForVerifying.size(); i++) {
        if (keepForVerifying[i]) {
            cumulativeSum++;
        }
        if (keepForPatterns[i]) {
            pattStars.push_back(cumulativeSum);
        }
    }

    // This is what's tricky about this step
    // We update starTable one last time, but also prepare pattStars for
    // database generation

    numEntries = (int)starTable.size();
    // size = 4358 = keptForPattCount
    // pattStars contains ints corresponding to star IDs

    std::cout << "For pattern matching with at most " << pattStarsPerFOV
              << " stars per FOV and no doubles: " << keepForPattCount
              << std::endl;
    std::cout << "For verification with at most " << verificationStarsPerFOV
              << " stars per FOV and no doubles: " << numEntries << std::endl
              << std::endl;
    std::cout << "Building temporary hash table for finding pattern neighbours"
              << std::endl;

    std::map<ShortVec3, std::set<short>> tempCoarseSkyMap;
    const int tempBins = 4;
    for (short starID : pattStars) {
        ShortVec3 hash = {
            (short)((starTable[starID].x + 1) * (float)tempBins),
            (short)((starTable[starID].y + 1) * (float)tempBins),
            (short)((starTable[starID].z + 1) * (float)tempBins),
        };
        tempCoarseSkyMap[hash].insert(starID);
    }

    // testing - looks good,  9/12/2022
    // for (auto const &pair : tempCoarseSkyMap) {
    //     std::cout << pair.first[0] << ", " << pair.first[1] << ", "
    //               << pair.first[2] << std::endl;
    //     for (auto ele : pair.second) {
    //         std::cout << ele << ", " << std::endl;
    //     }
    // }

    auto tempGetNearbyStars = [&tempCoarseSkyMap, &starTable](
                                  FloatVec3 vector, float radius, int ind) {
        std::vector<std::vector<int>> hcSpace;
        for (float x : vector) {
            std::vector<int> range;
            int lo = int((x + 1 - radius) * tempBins);
            lo = std::max(lo, 0);
            int hi = int((x + 1 + radius) * tempBins);
            hi = std::min(hi + 1, 2 * tempBins);
            range.push_back(lo);
            range.push_back(hi);
            hcSpace.push_back(range);
        }

        std::vector<short> nearbyStarIDs;

        for (int a = hcSpace[0][0]; a < hcSpace[0][1]; a++) {
            for (int b = hcSpace[1][0]; b < hcSpace[1][1]; b++) {
                for (int c = hcSpace[2][0]; c < hcSpace[2][1]; c++) {
                    ShortVec3 code{short(a), short(b), short(c)};
                    // possible that code is not in map?
                    if (tempCoarseSkyMap.count(code) == 0) {
                        continue;
                    }
                    for (short starID : tempCoarseSkyMap[code]) {
                        float dotProd = vector[0] * starTable[starID].x +
                                        vector[1] * starTable[starID].y +
                                        vector[2] * starTable[starID].z;

                        // if(ind == 1507 && starID == 3611){
                        //     nearbyStarIDs.push_back(3611);
                        // } // makes no difference
                        if (dotProd > std::cos(radius)) {
                            // if (ind == 1507) {
                            //     std::cout << "accepted, " << dotProd
                            //               << std::endl;
                            // }
                            // missing starID = 3611
                            // dotProd = cos(radius), NOT ACCEPTED 0.978148,
                            // 0.978148
                            nearbyStarIDs.push_back(starID);
                        }
                    }
                }
            }
        }

        return nearbyStarIDs;
    };

    typedef std::array<short, pattSize> Pattern;
    std::vector<Pattern> pattList;
    Pattern patt{0, 0, 0, 0};

    // std::cout << "patt stars size: " << pattStars.size() << std::endl;
    // 4358 stars in pattList

    for (int ind = 0; ind < (int)pattStars.size(); ind++) {
        // for(int ind = (int)pattStars.size()-1; ind >=0; ind--){
        short firstStarID = pattStars[ind];
        patt[0] = firstStarID;  // added

        // std::cout << "Pattern: " << ind << std::endl;
        // for (short x : patt) {
        //     std::cout << x << ", ";
        // }
        // std::cout << std::endl;
        // BUG: wrong here

        FloatVec3 firstStarVec = {
            starTable[firstStarID].x,
            starTable[firstStarID].y,
            starTable[firstStarID].z,
        };

        ShortVec3 hashCode = {
            short((firstStarVec[0] + 1) * (float)tempBins),
            short((firstStarVec[1] + 1) * (float)tempBins),
            short((firstStarVec[2] + 1) * (float)tempBins),
        };

        // std::cout << "hash code" << std::endl;
        // std::cout << hashCode[0] << ", " << hashCode[1] << ", " <<
        // hashCode[2] << std::endl; looks ok

        // tempCoarseSkyMap[hashCode].erase(
        //     tempCoarseSkyMap[hashCode].find(firstStarID));
        // Makes no difference?
        tempCoarseSkyMap[hashCode].erase(patt[0]);

        auto nearbyStars = tempGetNearbyStars(firstStarVec, maxFov, ind);
        // TODO: bug with tempGetNearbyStars, nearby stars is too low
        const int n = (int)nearbyStars.size();

        // std::cout << ind << ": " << n << std::endl;
        // BUG: for ind=1507, n = 26. However, in tetra3, n = 27
        // std::cout << "NEARBY STARS : " << n << std::endl;

        // BUG: doesn't even run pattern for patterns: 0, 1
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    patt[1] = nearbyStars[i];
                    patt[2] = nearbyStars[j];
                    patt[3] = nearbyStars[k];
                    // std::cout << "Pattern again " << patt[0] << ", " <<
                    // patt[1] << ", " << patt[2]
                    //           << ", " << patt[3] << std::endl;

                    // retrieve vectors of the stars in pattern
                    // check distance between every pair of stars in pattern
                    bool pattFits = true;
                    for (int pair1 = 0; pair1 < pattSize; pair1++) {
                        for (int pair2 = pair1 + 1; pair2 < pattSize; pair2++) {
                            float dotProd = starTable[patt[pair1]].x *
                                                starTable[patt[pair2]].x +
                                            starTable[patt[pair1]].y *
                                                starTable[patt[pair2]].y +
                                            starTable[patt[pair1]].z *
                                                starTable[patt[pair2]].z;
                            if (dotProd <= std::cos(maxFov)) {
                                pattFits = false;
                                break;
                            }
                        }
                    }
                    if (pattFits) {
                        // std::cout << patt[0] << ", " << patt[1] << ", "
                        //           << patt[2] << ", " << patt[3] <<
                        //           std::endl;
                        pattList.push_back(patt);
                    }
                }
            }
        }  // end of combination for loop
    }
    std::cout << "Found " << pattList.size() << " patterns" << std::endl;
    // BUG: missing 131 patterns

    int catalogLength = 2 * (int)pattList.size();
    std::vector<Pattern> pattCatalog(catalogLength);

    //
    for (Pattern patt : pattList) {
        std::vector<float> pattEdgeLengths;
        for (int i = 0; i < pattSize; i++) {
            StarTableEntry star1 = starTable[patt[i]];
            for (int j = i + 1; j < pattSize; j++) {
                // calculate distance between vectors
                StarTableEntry star2 = starTable[patt[j]];
                float edgeLen = std::sqrt(std::pow(star2.x - star1.x, 2) +
                                          std::pow(star2.y - star1.y, 2) +
                                          std::pow(star2.z - star1.z, 2));
                pattEdgeLengths.push_back(edgeLen);
            }
        }
        std::sort(pattEdgeLengths.begin(), pattEdgeLengths.end());

        float pattLargestEdge =
            pattEdgeLengths[(int)(pattEdgeLengths.size() - 1)];
        std::vector<float> pattEdgeRatios;
        for (int i = 0; i < (int)pattEdgeLengths.size() - 1;
             i++) {  // size()-1, since we ignore the largest edge
            pattEdgeRatios.push_back(pattEdgeLengths[i] / pattLargestEdge);
        }

        std::vector<int> key;
        for (float edgeRatio : pattEdgeRatios) {
            key.push_back(int(edgeRatio * pattBins));
        }

        int hashIndex = KeyToIndex(key, pattBins, catalogLength);
        long long offset = 0;
        while (true) {
            int index = int(hashIndex + std::pow(offset, 2)) % catalogLength;
            offset++;
            if (pattCatalog[index][0] == 0) {
                pattCatalog[index] = patt;
                break;
            }
        }
    }
    // DONE WITH EVERYTHING

    std::cout << "Catalog size: " << pattCatalog.size() << std::endl;
    std::cout << "Star table size: " << starTable.size() << std::endl;

    std::ofstream pattCatFile("pattCatalogOurs9-22.dat", std::ios::binary);
    // int pattInd = 0;
    for (Pattern patt : pattCatalog) {
        for(int i = 0; i < pattSize; i++){
            pattCatFile.write((char *)&patt[i], sizeof(short));
        }

        // std::cout << pattInd << ": " << patt[0] << ", " << patt[1] << ", "
        //           << patt[2] << ", " << patt[3] << std::endl;
        // pattInd++;
    }
    pattCatFile.close();

    std::ofstream starTableFile("starTableOurs9-22.dat", std::ios::binary);
    for(StarTableEntry ste: starTable){
        starTableFile.write((char*)&ste.ra, sizeof(float));
        starTableFile.write((char *)&ste.dec, sizeof(float));
        starTableFile.write((char *)&ste.x, sizeof(float));
        starTableFile.write((char *)&ste.y, sizeof(float));
        starTableFile.write((char *)&ste.z, sizeof(float));
        starTableFile.write((char *)&ste.magnitude, sizeof(float));
        starTableFile.write((char *)&ste.id, sizeof(float));
    }


    starTableFile.close();

}