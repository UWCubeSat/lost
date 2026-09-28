#include "databases-builder.hpp"

#include <bitset>
#include <iostream>
#include <string>

#include "attitude-utils.hpp"
#include "databases.hpp"
#include "decimal.hpp"
#include "io-util.hpp"
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

/// Build a star database from the catalog and write it to the output path.
void DatabaseBuild(const DatabaseOptions &values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(), (int) (values.minMag * 100), values.maxStars, DegToRad(values.minSeparation));
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseDescriptor dbEntries = GenerateDatabases(narrowedCatalog, values);
    SerializeContext ser = serFromDbValues(values);

    // Build the flags word. Currently the only flag indicates whether the
    // database was built with single-precision (float) or double-precision decimals.
    uint32_t dbFlags = 0;
    dbFlags |= typeid(decimal) == typeid(float) ? MULTI_DB_FLOAT_FLAG : 0;

    SerializeMultiDatabase(&ser, dbEntries, dbFlags);

    std::cerr << "Generated database with " << ser.buffer.size() << " bytes" << std::endl;
    std::cerr << "Database flagged with " << std::bitset<8*sizeof(dbFlags)>(dbFlags) << std::endl;

    UserSpecifiedOutputStream pos = UserSpecifiedOutputStream(values.outputPath, true);
    pos.Stream().write((char *) ser.buffer.data(), ser.buffer.size());
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
