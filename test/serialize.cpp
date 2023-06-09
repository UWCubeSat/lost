// Tests for serialization and multidatabases

#include <vector>

#include <catch.hpp>

#include "serialize-helpers.hpp"

using namespace lost; // NOLINT

TEST_CASE("Simple serialization, deserialization of primitives", "[fast] [serialize]") {
    int64_t val64 = 27837492938;
    float valFloat = 23.71728;
    SerializeContext ser;
    SerializePrimitive<int64_t>(&ser, val64);
    SerializePrimitive<float>(&ser, valFloat);
    DeserializeContext des(ser.buffer.data());
    int64_t deserializedVal64 = DeserializePrimitive<int64_t>(&des);
    float deserializedFloat = DeserializePrimitive<float>(&des);
    CHECK(val64 == deserializedVal64);
    CHECK(valFloat == deserializedFloat);
}

TEST_CASE("Endian-swapped serialization, deserialization of primitives", "[fast] [serialize]") {
    int64_t val64 = 27837492938;
    float valFloat = 23.71728;
    SerializeContext ser1(true, true);
    SerializePrimitive<int64_t>(&ser1, val64);
    SerializePrimitive<float>(&ser1, valFloat);
    DeserializeContext des(ser1.buffer.data());
    int64_t deserializedVal64 = DeserializePrimitive<int64_t>(&des);
    float deserializedValFloat = DeserializePrimitive<float>(&des);
    CHECK(val64 != deserializedVal64);
    CHECK(valFloat != deserializedValFloat);
    // but if we serialize it again, it should be back to normal!

    SerializeContext ser2(true, true);
    SerializePrimitive<int64_t>(&ser2, deserializedVal64);
    SerializePrimitive<float>(&ser2, deserializedValFloat);
    DeserializeContext des2(ser2.buffer.data());
    int64_t redeserializedVal64 = DeserializePrimitive<int64_t>(&des2);
    float redeserializedValFloat = DeserializePrimitive<float>(&des2);
    CHECK(val64 == redeserializedVal64);
    CHECK(valFloat == redeserializedValFloat);
}

TEST_CASE("Endian-swapped floats only", "[fast] [serialize]") {
    int64_t val64 = 27837492938;
    float valFloat = 23.71728;
    SerializeContext ser1(false, true);
    SerializePrimitive<int64_t>(&ser1, val64);
    SerializePrimitive<float>(&ser1, valFloat);
    DeserializeContext des(ser1.buffer.data());
    int64_t deserializedVal64 = DeserializePrimitive<int64_t>(&des);
    float deserializedValFloat = DeserializePrimitive<float>(&des);
    CHECK(val64 == deserializedVal64);
    CHECK(valFloat != deserializedValFloat);

    SerializeContext ser2(false, true);
    SerializePrimitive<float>(&ser2, deserializedValFloat);
    DeserializeContext des2(ser2.buffer.data());
    float redeserializedValFloat = DeserializePrimitive<float>(&des2);
    CHECK(valFloat == redeserializedValFloat);
}
