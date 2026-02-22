#ifndef DATABASES_BUILDER_H
#define DATABASES_BUILDER_H

#include <string>
#include <iostream>

#include "databases.hpp"
#include "camera.hpp"
#include "star-utils.hpp"
#include "decimal.hpp"

namespace lost {

/// Commannd line options when using the `database` command.
class DatabaseOptions {
public:
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
    type prop = defaultVal;
#include "database-options.hpp"
#undef LOST_CLI_OPTION
};

SerializeContext serFromDbValues(const DatabaseOptions &values);

/// Appropriately create descriptors for all requested databases according to command-line options.
/// @sa SerializeMultiDatabase
MultiDatabaseDescriptor GenerateDatabases(const Catalog &, const DatabaseOptions &values);

/// Build a star database from the catalog and write it to the output path.
void DatabaseBuild(const DatabaseOptions &values);

std::ostream &operator<<(std::ostream &, const Camera &);

}

#endif
