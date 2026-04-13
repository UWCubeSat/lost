#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <utility>

#ifdef LOST_USE_ETL_CONTAINERS
#include "etl/memory.h"
#include "etl/pool.h"
#include "etl/vector.h"
#else
#include <vector>
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
