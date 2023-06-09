#include "utils.hpp"

#include <algorithm>
#include <vector>

#include "databases.hpp"

namespace lost {

bool AreStarIdentifiersEquivalent(const StarIdentifiers &ids1, const StarIdentifiers &ids2) {
    if (ids1.size() != ids2.size()) {
        return false;
    }

    for (const StarIdentifier &id1 : ids1) {
        if (std::find(ids2.cbegin(), ids2.cend(), id1) == ids2.cend()) {
            return false;
        }
    }
    return true;
}

}
