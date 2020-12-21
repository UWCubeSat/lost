#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <vector>

#include "catalog-generic.hpp"

class DatabaseBuilder {
public:
    virtual unsigned char *Go(int *length) const;
};

/**
 For a K-vector database with `n' stars and `m' distance bins
     | size          | description         |
     |---------------+---------------------|
     | sizeof(int)   | stores `n'          |
     | n*sizeof(int) | stores star indices, sorted  |
     |               |                     |
 */

class KVectorDatabaseBuilder {
public:
    KVectorDatabaseBuilder(long minDistance, long maxDistance)
        : minDistance(minDistance), maxDistance(maxDistance) { };
    unsigned char *Go(int *length) const;
private:
    long minDistance; // millionth of a degree-
    long maxDistance;
};

KVectorDatabaseBuilder PromptKVectorDatabaseBuilder();

class KVectorDatabase {
    std::vector<int> StarsWithDistance(long minDistance, long maxDistance);
private:
    unsigned char *bytes;
};

#endif
