#ifndef CATALOG_GENERIC_H
#define CATALOG_GENERIC_H

namespace lost {

class CatalogStar {
public:
    CatalogStar(long raj2000, long dej2000, int magnitude, bool weird) :
        raj2000(raj2000), dej2000(dej2000), magnitude(magnitude), weird(weird) { }
    long raj2000;           // *10^-6, right ascension
    long dej2000;           // *10^-6, declination
    int  magnitude;         // *10^-2
    bool weird;             // nonzero for binary, etc
};

}

#endif
