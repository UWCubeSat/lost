#ifndef CONTAINERS_H
#define CONTAINERS_H

// ETL Configuration for embedded systems
// Maximum entries for the narrowed star catalog.
// Matches the default value of the --max-stars database option.
#ifndef LOST_ETL_MAX_CATALOG_STARS
#define LOST_ETL_MAX_CATALOG_STARS 10000
#endif

// Maximum pixel count for generated/runtime images in ETL mode.
// Matches the default generated image size: 1024 * 1024.
#ifndef LOST_ETL_MAX_IMAGE_PIXELS
#define LOST_ETL_MAX_IMAGE_PIXELS 1048576
#endif

// Maximum number of centroid/star records kept in memory for one frame.
#ifndef LOST_ETL_MAX_STARS
#define LOST_ETL_MAX_STARS 4096
#endif

// Maximum temporary serialization buffer size in bytes.
#ifndef LOST_ETL_MAX_SERIALIZE_BUFFER_BYTES
#define LOST_ETL_MAX_SERIALIZE_BUFFER_BYTES 4194304
#endif

// Capacity for command pipeline inputs (single image by default).
#ifndef LOST_ETL_MAX_PIPELINE_INPUTS
#define LOST_ETL_MAX_PIPELINE_INPUTS 4
#endif

// Pool sizes for pointer-heavy runtime objects used in ETL mode.
#ifndef LOST_ETL_MAX_CENTROID_ALGO_OBJECTS
#define LOST_ETL_MAX_CENTROID_ALGO_OBJECTS 2
#endif

#ifndef LOST_ETL_MAX_STAR_ID_ALGO_OBJECTS
#define LOST_ETL_MAX_STAR_ID_ALGO_OBJECTS 2
#endif

#ifndef LOST_ETL_MAX_ATTITUDE_ALGO_OBJECTS
#define LOST_ETL_MAX_ATTITUDE_ALGO_OBJECTS 2
#endif

#ifndef LOST_ETL_MAX_PIPELINE_STARS_OBJECTS
#define LOST_ETL_MAX_PIPELINE_STARS_OBJECTS 8
#endif

#ifndef LOST_ETL_MAX_PIPELINE_STAR_IDS_OBJECTS
#define LOST_ETL_MAX_PIPELINE_STAR_IDS_OBJECTS 8
#endif

#ifndef LOST_ETL_MAX_PIPELINE_ATTITUDE_OBJECTS
#define LOST_ETL_MAX_PIPELINE_ATTITUDE_OBJECTS 8
#endif

#ifndef LOST_ETL_MAX_IDENTIFIED_STARS_IN_RANGE
#define LOST_ETL_MAX_IDENTIFIED_STARS_IN_RANGE 64
#endif

#ifndef LOST_ETL_MAX_MULTI_DATABASE_ENTRIES
#define LOST_ETL_MAX_MULTI_DATABASE_ENTRIES 8
#endif

#ifndef LOST_ETL_MAX_KVECTOR_DISTANCE_BINS_PLUS_ONE
#define LOST_ETL_MAX_KVECTOR_DISTANCE_BINS_PLUS_ONE 10001
#endif

#ifndef LOST_ETL_MAX_PAIR_DISTANCE_PAIRS
#define LOST_ETL_MAX_PAIR_DISTANCE_PAIRS 500000
#endif

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>

#ifdef LOST_USE_ETL_CONTAINERS
#include "etl/memory.h"
#include "etl/pool.h"
#include "etl/vector.h"
#include "etl/bitset.h"
#else
#include <vector>
#include <bitset>
#endif
namespace lost {

#ifdef LOST_USE_ETL_CONTAINERS
struct pool_deleter {
    void *pool = nullptr;
    void *object = nullptr;
    void (*destroy_fn)(void *, void *) = nullptr;

    template <typename PointerType>
    void operator()(PointerType *ptr) const {
        (void)ptr;
        if ((destroy_fn != nullptr) && (object != nullptr)) {
            destroy_fn(pool, object);
        }
    }

    template <typename PoolType, typename ObjectType>
    static pool_deleter Make(PoolType &pool_ref, ObjectType *object_ptr) {
        return {&pool_ref, object_ptr, [](void *pool_ptr, void *stored_object) {
                    static_cast<PoolType *>(pool_ptr)->destroy(static_cast<ObjectType *>(stored_object));
                }};
    }

};

template <typename T, size_t N>
using vector = etl::vector<T, N>;

template <size_t N>
using bitset = etl::bitset<N>;

template <typename T, size_t N>
using pool = etl::pool<T, N>;

template <typename T, size_t N = 0>
using unique_ptr = etl::unique_ptr<T, pool_deleter>;

template <typename T, size_t N = 0, typename PoolType, typename... Args>
unique_ptr<T, N> make_unique(PoolType &pool_ref, Args &&...args) {
    T *object = pool_ref.create(std::forward<Args>(args)...);
    if (object == nullptr) {
        std::cerr << "ERROR: ETL pool exhausted while allocating object." << std::endl;
        std::exit(1);
    }
    return unique_ptr<T, N>(object, pool_deleter::Make(pool_ref, object));
}

template <typename BaseType, size_t N = 0, typename DerivedType, typename PoolType, typename... Args>
unique_ptr<BaseType, N> make_unique_base(PoolType &pool_ref, Args &&...args) {
    DerivedType *object = pool_ref.create(std::forward<Args>(args)...);
    if (object == nullptr) {
        std::cerr << "ERROR: ETL pool exhausted while allocating polymorphic object." << std::endl;
        std::exit(1);
    }
    return unique_ptr<BaseType, N>(object, pool_deleter::Make(pool_ref, object));
}

#else
struct dummy_pool {
};

template <typename T, size_t N>
using vector = std::vector<T>;

template <size_t N>
using bitset = std::bitset<N>;

template <typename T, size_t N>
using pool = dummy_pool;

template <typename T, size_t N = 0>
using unique_ptr = std::unique_ptr<T>;

template <typename T, size_t N = 0, typename PoolType, typename... Args>
unique_ptr<T, N> make_unique(PoolType &pool_ref, Args &&...args) {
    (void)pool_ref;
    return std::make_unique<T>(std::forward<Args>(args)...);
}

template <typename BaseType, size_t N = 0, typename DerivedType, typename PoolType, typename... Args>
unique_ptr<BaseType, N> make_unique_base(PoolType &pool_ref, Args &&...args) {
    (void)pool_ref;
    std::unique_ptr<DerivedType> derived = std::make_unique<DerivedType>(std::forward<Args>(args)...);
    return unique_ptr<BaseType, N>(std::move(derived));
}

#endif

inline void EtlRuntimeBoundCheck(std::size_t value, std::size_t maxAllowed, const char *label) {
#ifdef LOST_USE_ETL_CONTAINERS
    if (value > maxAllowed) {
        std::cerr << "ERROR: " << label << " exceeded ETL bound. value=" << value
                  << " max=" << maxAllowed << std::endl;
        std::exit(1);
    }
#else
    (void)value;
    (void)maxAllowed;
    (void)label;
#endif
}

}

#endif
