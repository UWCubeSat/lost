/**
 * Helpers to serialize and deserialize arbitrary data types to disk
 *
 * The serialization and deserialization helpers here assume that (a) integers and floating point
 * numbers are stored in the same format on the source and target systems except for endianness, and
 * (b) the target system requires no greater than n-byte alignment for n-byte values (eg, an int32_t
 * only needs 4-byte alignment, not 8-byte), and (c) the database itself will be aligned at a
 * multiple of the largest size of any value stored within (because we only align values relative to
 * the start of the database, so everything breaks if the database itself is misaligned).
 *
 * Generally, the "bulk" of any database will not be explicitly deserialized into a data structure
 * in memory, but instead will just effectively memory-mapped, by just storing a pointer to a
 * certain offset into the database. The serialization functions in this file can still be used to
 * create such a "bulk" section, but deserialization must be handled manually.
 */

#ifndef SERIALIZE_HELPERS_H
#define SERIALIZE_HELPERS_H

#include <utility>
#include <vector>
#include <string.h>

#ifndef LOST_DATABASE_SOURCE_INTEGER_ENDIANNESS
#define LOST_DATABASE_SOURCE_INTEGER_ENDIANNESS lost::Endianness::LittleEndian
#endif

#ifndef LOST_DATABASE_SOURCE_FLOAT_ENDIANNESS
#define LOST_DATABASE_SOURCE_FLOAT_ENDIANNESS lost::Endianness::LittleEndian
#endif

#ifndef LOST_DATABASE_TARGET_INTEGER_ENDIANNESS
#define LOST_DATABASE_TARGET_INTEGER_ENDIANNESS lost::Endianness::LittleEndian
#endif

#ifndef LOST_DATABASE_TARGET_FLOAT_ENDIANNESS
#define LOST_DATABASE_TARGET_FLOAT_ENDIANNESS lost::Endianness::LittleEndian
#endif

namespace lost {

enum class Endianness {
    LittleEndian,
    BigEndian,
};

class DeserializeContext {
public:
    DeserializeContext(const unsigned char *buffer) : buffer(buffer), cursor(buffer) { };

    size_t GetOffset() const {
        return cursor - buffer;
    }

    void MoveForward(size_t howMuch) {
        cursor += howMuch;
    }

    const unsigned char *GetCursor() {
        return cursor;
    }

private:
    const unsigned char *buffer; /// the start of the buffer
    const unsigned char *cursor; /// the current location of the "read head" into the buffer
};

/// Unconditionally swap the endianness of a value (uses sizeof T).
template <typename T>
void SwapEndianness(T *buffer) {
    char *charBuffer = (char *)buffer;
    for (int i = 0; i < (int)(sizeof(T)/2); i++) {
        std::swap(charBuffer[i], charBuffer[sizeof(T)-1-i]);
    }
}

/// Swap the endianness of a value if necessary. Uses
/// LOST_DATABASE_{SOURCE,TARGET}_INTEGER_ENDIANNESS to determine to switch all values but float and
/// double, which use LOST_DATABASE_{SOURCE,TARGET}_FLOAT_ENDIANNESS.
template <typename T>
void SwapEndiannessIfNecessary(T *buffer) {
    if (LOST_DATABASE_SOURCE_INTEGER_ENDIANNESS != LOST_DATABASE_TARGET_INTEGER_ENDIANNESS) {
        SwapEndianness<T>(buffer);
    }
}

// template specializations

template <>
inline void SwapEndiannessIfNecessary<float>(float *buffer) {
    if (LOST_DATABASE_SOURCE_FLOAT_ENDIANNESS != LOST_DATABASE_TARGET_FLOAT_ENDIANNESS) {
        SwapEndianness<float>(buffer);
    }
}

template <>
inline void SwapEndiannessIfNecessary<double>(double *buffer) {
    if (LOST_DATABASE_SOURCE_FLOAT_ENDIANNESS != LOST_DATABASE_TARGET_FLOAT_ENDIANNESS) {
        SwapEndianness<double>(buffer);
    }
}

/// Move the cursor forward past any padding that would appear before a value of type T
template <typename T>
void DeserializePadding(DeserializeContext *des) {
    des->MoveForward((sizeof(T) - des->GetOffset()%sizeof(T))%sizeof(T));
}

template <typename T>
T DeserializePrimitive(DeserializeContext *des) {
    DeserializePadding<T>(des);
    const T *result = (T *)des->GetCursor(); // endianness should have been taken care of during serialization
    des->MoveForward(sizeof(T));
    return *result;
}

/// return an array of items as a pointer. Will point into the buffer (mmap style).
template <typename T>
const T *DeserializeArray(DeserializeContext *des, long arrLength) {
    DeserializePadding<T>(des); // Perhaps we should always 8-align arrays? Does that possibly make
                                // any offset calculation or accesses faster?
    const T *result = (T *)des->GetCursor();
    des->MoveForward(sizeof(T)*arrLength);
    return result;
}

template <typename T>
void SerializePadding(std::vector<unsigned char> *vec) {
    for (int i = vec->size(); i%sizeof(T) != 0; i++) {
        vec->push_back(0);
    }
}

template <typename T>
void SerializePrimitive(std::vector<unsigned char> *vec, const T &val) {
    T tBuf;
    memcpy(&tBuf, &val, sizeof(T));
    SwapEndiannessIfNecessary<T>(&tBuf);
    unsigned char *buf = (unsigned char *)(&tBuf);
    SerializePadding<T>(vec);
    std::copy(buf, buf+sizeof(T), std::back_inserter(*vec));
}

}

#endif
