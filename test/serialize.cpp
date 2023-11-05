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

TEST_CASE("Padding", "[fast] [serialize]") {
    int8_t val8 = 23;
    int32_t val32 = 1234567;
    SerializeContext ser;
    SerializePrimitive<int8_t>(&ser, val8);
    SerializePrimitive<int32_t>(&ser, val32);
    CHECK(ser.buffer.size() == 8);
    CHECK(ser.buffer[0] == 23);
    CHECK(ser.buffer[1] == 0);
    CHECK(ser.buffer[2] == 0);
    CHECK(ser.buffer[3] == 0);

    DeserializeContext des(ser.buffer.data()+4);
    int32_t deserializedVal32 = DeserializePrimitive<int32_t>(&des);
    CHECK(val32 == deserializedVal32);
}

TEST_CASE("Array", "[fast] [serialize]") {
    SerializeContext ser;
    // serialize a single byte first, to ensure that array adds padding.
    SerializePrimitive<int8_t>(&ser, 42);
    for (int16_t i = 0; i < 16; i++) {
        SerializePrimitive<int16_t>(&ser, i);
    }

    DeserializeContext des(ser.buffer.data());
    int8_t firstByte = DeserializePrimitive<int8_t>(&des);
    CHECK(firstByte == 42);
    const int16_t *arr = DeserializeArray<int16_t>(&des, 8);
    CHECK(arr[0] == 0);
    CHECK(arr[1] == 1);
    int8_t ninthByte = DeserializePrimitive<int16_t>(&des);
    CHECK(ninthByte == 8);
}
