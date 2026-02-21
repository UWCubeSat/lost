#include "databases-builder.hpp"

#include <iostream>
#include <string>

#include "attitude-utils.hpp"
#include "databases.hpp"
#include "decimal.hpp"
#include "star-utils.hpp"

namespace lost {

SerializeContext serFromDbValues(const DatabaseOptions &values) {
    return SerializeContext(values.swapIntegerEndianness, values.swapDecimalEndianness);
}

MultiDatabaseDescriptor GenerateDatabases(const Catalog &catalog, const DatabaseOptions &values) {
    MultiDatabaseDescriptor dbEntries;

    SerializeContext catalogSer = serFromDbValues(values);
    // TODO decide why we have this inclMagnitude and inclName and if we should change that
    SerializeCatalog(&catalogSer, catalog, false, true);
    dbEntries.emplace_back(kCatalogMagicValue, catalogSer.buffer);

    if (values.kvector) {
        decimal minDistance = DegToRad(values.kvectorMinDistance);
        decimal maxDistance = DegToRad(values.kvectorMaxDistance);
        long numBins = values.kvectorNumDistanceBins;
        SerializeContext ser = serFromDbValues(values);
        SerializePairDistanceKVector(&ser, catalog, minDistance, maxDistance, numBins);
        dbEntries.emplace_back(PairDistanceKVectorDatabase::kMagicValue, ser.buffer);
    } else {
        std::cerr << "No database builder selected -- no database generated." << std::endl;
        exit(1);
    }

    return dbEntries;
}

/// Print information about the camera in machine and human-readable form.
std::ostream &operator<<(std::ostream &os, const Camera &camera) {
    os << "camera_focal_length " << camera.FocalLength() << std::endl
       << "camera_fov " << camera.Fov() << std::endl
       << "camera_resolution_x " << camera.XResolution() << std::endl
       << "camera_resolution_y " << camera.YResolution() << std::endl;
    // TODO: principal point
    return os;
}

} // namespace lost
