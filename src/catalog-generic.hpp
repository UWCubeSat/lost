#ifndef CATALOG_GENERIC_H
#define CATALOG_GENERIC_H

#include <vector>
#include <string>

namespace lost {

class CatalogStar {
public:
    CatalogStar(float raj2000, float dej2000, int magnitude, bool weird, std::string name) :
        raj2000(raj2000), dej2000(dej2000), magnitude(magnitude), weird(weird), name(name) { }
    float raj2000;           // *10^-6, right ascension
    float dej2000;           // *10^-6, declination
    int  magnitude;         // *10^-2
    bool weird;             // nonzero for binary, etc
    std::string name;
};

typedef std::vector<CatalogStar> Catalog;

}

#endif
