#ifndef ETL_CONFIG_H
#define ETL_CONFIG_H

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

// Maximum number of star identifiers produced for one frame.
#ifndef LOST_ETL_MAX_STAR_IDENTIFIERS
#define LOST_ETL_MAX_STAR_IDENTIFIERS LOST_ETL_MAX_STARS
#endif

// Maximum number of pair-distance query candidates used during star-id.
#ifndef LOST_ETL_MAX_PAIR_QUERY_RESULTS
#define LOST_ETL_MAX_PAIR_QUERY_RESULTS 16384
#endif

// Maximum temporary serialization buffer size in bytes.
#ifndef LOST_ETL_MAX_SERIALIZE_BUFFER_BYTES
#define LOST_ETL_MAX_SERIALIZE_BUFFER_BYTES 4194304
#endif

// Capacity for command pipeline inputs (single image by default).
#ifndef LOST_ETL_MAX_PIPELINE_INPUTS
#define LOST_ETL_MAX_PIPELINE_INPUTS 4
#endif

// Capacity for command pipeline outputs.
#ifndef LOST_ETL_MAX_PIPELINE_OUTPUTS
#define LOST_ETL_MAX_PIPELINE_OUTPUTS LOST_ETL_MAX_PIPELINE_INPUTS
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

// Temporary container capacities for image generation and comparison paths.
#ifndef LOST_ETL_MAX_GENERATED_STARS
#define LOST_ETL_MAX_GENERATED_STARS LOST_ETL_MAX_CATALOG_STARS
#endif

#ifndef LOST_ETL_MAX_CLOSEST_STAR_CANDIDATES
#define LOST_ETL_MAX_CLOSEST_STAR_CANDIDATES LOST_ETL_MAX_STARS
#endif

#ifndef LOST_ETL_MAX_CENTROID_COMPARISONS
#define LOST_ETL_MAX_CENTROID_COMPARISONS LOST_ETL_MAX_PIPELINE_OUTPUTS
#endif

#ifndef LOST_ETL_MAX_THIRD_STAR_CANDIDATES
#define LOST_ETL_MAX_THIRD_STAR_CANDIDATES LOST_ETL_MAX_PAIR_QUERY_RESULTS
#endif

#ifndef LOST_ETL_MAX_UNIDENTIFIED_CENTROIDS
#define LOST_ETL_MAX_UNIDENTIFIED_CENTROIDS LOST_ETL_MAX_STARS
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

#endif
