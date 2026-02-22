#ifndef IO_UTIL_H
#define IO_UTIL_H

#include <string>
#include <iostream>

#include "star-utils.hpp"

namespace lost {

const char kNoDefaultArgument = 0;

/// An output stream which might be a file or stdout
class UserSpecifiedOutputStream {
public:
    explicit UserSpecifiedOutputStream(std::string filePath, bool isBinary);
    ~UserSpecifiedOutputStream();

    /// return the inner output stream, suitable for use with <<
    std::ostream &Stream() { return *stream; };

private:
    bool isFstream;
    std::ostream *stream;
};

// use the environment variable LOST_BSC_PATH, or read from ./bright-star-catalog.tsv
const Catalog &CatalogRead();

}

#endif
