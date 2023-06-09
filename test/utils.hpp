#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "databases.hpp"

namespace lost {

/// simple O(n^2) check
bool AreStarIdentifiersEquivalent(const StarIdentifiers &, const StarIdentifiers &);

}

#endif
