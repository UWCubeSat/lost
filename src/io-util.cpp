#include "io-util.hpp"

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

#include "attitude-utils.hpp"
#include "decimal.hpp"
#include "star-utils.hpp"

namespace lost {

/// Create a PromptedOutputStream which will output to the given file.
UserSpecifiedOutputStream::UserSpecifiedOutputStream(std::string filePath, bool isBinary) {
    if (isBinary && isatty(fileno(stdout)) && (filePath == "stdout" || filePath == "-")) {
        std::cerr << "WARNING: output contains binary contents. Not printed to terminal." << std::endl;
        filePath = "/dev/null";
    }

    if (filePath == "stdout" || filePath == "-") {
        stream = &std::cout;
        isFstream = false;
    } else {
        std::fstream *fs = new std::fstream();
        fs->open(filePath, std::fstream::out);
        stream = fs;
        isFstream = true;
    }
}

UserSpecifiedOutputStream::~UserSpecifiedOutputStream() {
    if (isFstream) {
        delete stream;
    }
}

/// Parse the bright star catalog from the TSV file on disk.
std::vector<CatalogStar> BscParse(std::string tsvPath) {
    std::vector<CatalogStar> result;
    FILE *file;
    decimal raj2000, dej2000;
    int magnitudeHigh, magnitudeLow, name;
    char weird;

    file = fopen(tsvPath.c_str(), "r");

    if (file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        exit(1); // TODO: do we want any other error handling?
        return result;
    }

    #ifdef LOST_FLOAT_MODE
        std::string format = "%f|%f|%d|%c|%d.%d";
    #else
        std::string format = "%lf|%lf|%d|%c|%d.%d";
    #endif

    while (EOF != fscanf(file, format.c_str(),
                         &raj2000, &dej2000,
                         &name, &weird,
                         &magnitudeHigh, &magnitudeLow)) {
        result.push_back(CatalogStar(DegToRad(raj2000),
                                     DegToRad(dej2000),
                                     magnitudeHigh*100 + (magnitudeHigh < 0 ? -magnitudeLow : magnitudeLow),
                                     name));
    }

    fclose(file);
    assert(result.size() > 9000); // basic sanity check
    return result;
}

#ifndef DEFAULT_BSC_PATH
#define DEFAULT_BSC_PATH "bright-star-catalog.tsv"
#endif

/// Read and parse the full catalog from disk. If called multiple times, will re-use the first result.
const Catalog &CatalogRead() {
    static bool readYet = false;
    static std::vector<CatalogStar> catalog;

    if (!readYet) {
        readYet = true;
        char *tsvPath = getenv("LOST_BSC_PATH");
        catalog = BscParse(tsvPath ? tsvPath : DEFAULT_BSC_PATH);
        // perform essential narrowing
        // remove all stars with exactly the same position as another, keeping the one with brighter magnitude
        std::sort(catalog.begin(), catalog.end(), [](const CatalogStar &a, const CatalogStar &b) {
            return a.spatial.x < b.spatial.x;
        });
        for (int i = catalog.size()-1; i > 0; i--) {
            if ((catalog[i].spatial - catalog[i-1].spatial).Magnitude() < DECIMAL(5e-5)) { // 70 stars removed at this threshold.
                if (catalog[i].magnitude > catalog[i-1].magnitude) {
                    catalog.erase(catalog.begin() + i);
                } else {
                    catalog.erase(catalog.begin() + i - 1);
                }
            }
        }
    }
    return catalog;
}

} // namespace lost
