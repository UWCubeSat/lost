# LOST Architecture & Porting Strategy

## 1. Original C++ Source → Python Module Mapping

```
C++ File                          → Python Module          Purpose
─────────────────────────────────────────────────────────────────────────────
main.cpp                          → (covered in)           CLI entry point
                                    pipeline.py

pipeline-options.hpp              → config.py              All pipeline configuration
database-options.hpp              → config.py              All database configuration

io.hpp / io.cpp                   → pipeline.py            Pipeline orchestration
  - Pipeline class                → pipeline.py              Pipeline runner
  - PipelineInput hierarchy       → pipeline.py              Input abstraction
  - PipelineOutput                → pipeline.py              Output data
  - PipelineComparison            → validation.py            Compare & report
  - SetPipeline()                 → pipeline.py              Factory from config
  - CatalogRead()                 → database.py              BSC loading
  - FocalLengthFromOptions()      → camera.py                Focal length from config
  - CentroidComparison            → validation.py            Centroid comparison
  - StarIdsCompare                → validation.py            Star ID comparison
  - SurfacePlot                   → (skip, Cairo GUI)        Image annotation

camera.hpp / camera.cpp           → camera.py               Camera model
  - Camera::SpatialToCamera       → Camera.spatial_to_camera  3D→2D projection
  - Camera::CameraToSpatial       → Camera.camera_to_spatial  2D→3D unprojection
  - Camera::InSensor              → Camera.in_sensor          Pixel bounds check
  - FovToFocalLength              → camera.fov_to_focal_length
  - FocalLengthToFov              → camera.focal_length_to_fov

centroiders.hpp / centroiders.cpp → centroiding.py          Centroid detection
  - CentroidAlgorithm (base)      → CentroidAlgorithm          Abstract base
  - DummyCentroidAlgorithm        → DummyCentroider            Debug centroids
  - CenterOfGravityAlgorithm      → CogCentroider              CoG method
  - IterativeWeightedCoGAlgorithm → IterativeWeightedCogCentroider  IWCoG method
  - BasicThreshold                → centroiding.basic_threshold  Image threshold
  - BasicThresholdOnePass         → centroiding.basic_threshold_onepass
  - OtsusThreshold                → centroiding.otsus_threshold
  - CogHelper                     → (inlined in CogCentroider)
  - IWCoGHelper                   → (inlined in IWCoG)

star-utils.hpp / star-utils.cpp   → star_id.py (types)      + database.py
  - CatalogStar                   → CatalogStar               Catalog entry
  - Star                          → Star                      Centroid
  - StarIdentifier                → StarIdentifier            Star match pair
  - Catalog                       → list[CatalogStar]         Type alias
  - Stars                         → list[Star]                Type alias
  - StarIdentifiers               → list[StarIdentifier]      Type alias
  - NarrowCatalog                 → database.narrow_catalog   Catalog filtering
  - MagToBrightness               → star_id.mag_to_brightness
  - SerializeCatalog              → (skip, binary format)     DB serialization
  - DeserializeCatalog            → (skip, binary format)     DB deserialization

star-id.hpp / star-id.cpp         → star_id.py               Star identification
  - StarIdAlgorithm (base)        → StarIdAlgorithm            Abstract base
  - DummyStarIdAlgorithm          → DummyStarId                Debug star ID
  - GeometricVotingStarIdAlgorithm→ GeometricVotingStarId      GV algorithm
  - PyramidStarIdAlgorithm        → PyramidStarId              Pyramid algorithm
  - PairDistanceInvolvingIterator → star_id._PairDistanceInvolving
  - PairDistanceQueryToMap        → star_id._pairs_to_map
  - IdentifyThirdStar             → star_id._identify_third_star
  - IdentifyRemainingStarsPairDist→ star_id._identify_remaining_stars
  - IRUnidentifiedCentroid        → star_id._IRUnidentifiedCentroid
  - SelectNextUnidentifiedCentroid→ star_id._select_next_unidentified
  - AddToAllUnidentifiedCentroids → star_id._add_to_all_unidentified

star-id-private.hpp               → star_id.py               Internal helpers

attitude-utils.hpp / .cpp         → math_utils.py            Math primitives
  - Vec2                          → Vec2                       2D vector
  - Vec3                          → Vec3                       3D vector
  - Mat3                          → Mat3                       3×3 matrix
  - Quaternion                    → Quaternion                 Rotation quaternion
  - Attitude                      → Attitude                   Combined attitude
  - EulerAngles                   → EulerAngles                RA/Dec/Roll
  - SphericalToSpatial            → math_utils.spherical_to_spatial
  - SpatialToSpherical            → math_utils.spatial_to_spherical
  - SphericalToQuaternion         → math_utils.spherical_to_quaternion
  - QuaternionToDCM               → math_utils.quaternion_to_dcm
  - DCMToQuaternion               → math_utils.dcm_to_quaternion
  - Angle                         → math_utils.angle
  - AngleUnit                     → math_utils.angle_unit
  - RadToDeg / DegToRad           → math_utils.rad_to_deg / deg_to_rad
  - RadToArcSec / ArcSecToRad     → math_utils.rad_to_arcsec / arcsec_to_rad
  - DecimalModulo                 → math_utils.decimal_modulo

attitude-estimators.hpp / .cpp    → attitude.py              Attitude estimation
  - AttitudeEstimationAlgorithm   → AttitudeEstimator          Abstract base
  - DavenportQAlgorithm           → DavenportQAttitude         DQM method
  - TriadAlgorithm                → TriadAttitude              TRIAD method
  - QuestAlgorithm                → QuestAttitude              QUEST method
  - QuestCharPoly                 → attitude._quest_char_poly
  - QuestCharPolyPrime            → attitude._quest_char_poly_prime
  - QuestEigenvalueEstimator      → attitude._quest_eigenvalue
  - TriadCoordinateFrame          → attitude._triad_frame

databases.hpp / databases.cpp     → database.py              Database structures
  - KVectorIndex                  → KVectorIndex               K-vector index
  - KVectorIndex::QueryLiberal    → kvector.query_liberal
  - KVectorIndex::BinFor          → kvector._bin_for
  - PairDistanceKVectorDatabase   → PairDistanceKVectorDB      Pair-distance DB
  - FindPairsLiberal              → find_pairs_liberal
  - FindPairsExact                → find_pairs_exact
  - MultiDatabase                 → MultiDatabase              Multi-DB container
  - SerializeKVectorIndex         → (skip, binary)             DB serialization
  - SerializePairDistanceKVector  → (skip, binary)             DB serialization
  - CatalogToPairDistances        → _catalog_to_pair_distances

serialize-helpers.hpp             → database.py              Serialization helpers
  - SerializeContext              → (skip)                     Binary serialization
  - DeserializeContext            → (skip)                     Binary deserialization

decimal.hpp                       → config.py                 Decimal type config
                                                              (use float directly)

------------------ NOT PORTED (image generation / benchmarking) ----------------

io.cpp (GeneratedPipelineInput)   → (skip)                   Image generation
io.cpp (RandomAttitude)           → (skip)                   Random attitude
io.cpp (SurfacePlot)              → (skip)                   Cairo plotting
io.cpp (PipelineComparison plot)  → (skip)                   Cairo-based plots
io.cpp (BscParse)                 → database.py              TSV catalog parser
```

## 2. Algorithm Descriptions

### 2.1 Centroiding Algorithms

**Algorithm: BasicThreshold**
- Source: `centroiders.cpp:81-94`
- Purpose: Compute pixel intensity threshold to separate stars from background
- Math: `threshold = mean(pixels) + 5 * stddev(pixels)`
- Input: Grayscale image (uint8 array)
- Output: Single threshold value

**Algorithm: Center of Gravity (CoG)**
- Source: `centroiders.cpp:158-196`
- Purpose: Detect star centroids by flood-fill + weighted average
- Method:
  1. Compute threshold via BasicThreshold
  2. Scan pixels; when pixel ≥ threshold and unvisited, flood-fill connected region
  3. Reject regions touching image edge
  4. Compute centroid as intensity-weighted average of pixel coordinates
  5. Compute radius from bounding box
- Formulas:
  - `x_center = Σ(x_i * I_i) / Σ(I_i)`
  - `y_center = Σ(y_i * I_i) / Σ(I_i)`
- Complexity: O(width × height × region_area)
- Note: Adds 0.5 pixel offset to coordinates (pixel-center convention)

**Algorithm: Iterative Weighted CoG (IWCoG)**
- Source: `centroiders.cpp:247-325`
- Purpose: More precise centroid via Gaussian-weighted iteration
- Method:
  1. Same flood-fill as CoG
  2. Find pixel with max intensity as initial guess
  3. Compute FWHM, convert to Gaussian standard deviation
  4. Iteratively re-estimate centroid using Gaussian-weighted pixel values
  5. Stop when change < 0.0002 or after 100000 iterations
- Gaussian weight: `w = I_max * exp(-((x-x_g)^2 + (y-y_g)^2) / (2*σ^2))`
- Complexity: O(width × height × iterations)
- Note: LOST comment says "doesn't perform much better than CoG"

### 2.2 Star Identification Algorithms

**Algorithm: Geometric Voting**
- Source: `star-id.cpp:29-152`
- Purpose: Identify stars by voting for catalog matches
- Method:
  1. For each star i in image:
     a. For each other star j:
        - Compute great-circle distance θ_ij
        - Query pair-distance DB for catalog pairs within tolerance
        - Vote for each catalog star appearing in matching pairs
     b. Assign the catalog star with the most votes
  2. Verification phase:
     - For each pair of identified stars (i,j):
       - Compute catalog distance and image distance
       - If |catalog_dist - image_dist| < tolerance, verify both
     - Keep only stars with votes ≥ 75% of max verification votes
- Complexity: O(N² × M) where N=image stars, M=catalog pairs per query
- Failure modes: Few stars, high false-star rate, large tolerance

**Algorithm: Pyramid**
- Source: `star-id.cpp:571-770`
- Purpose: Robust 4-star pattern matching with probabilistic validation
- Method:
  1. Enumerate 4-star tuples (i,j,k,r) using interleaved index patterns
  2. For each tuple:
     a. Compute 6 inter-star distances (ij, ik, ir, jk, jr, kr)
     b. Compute 3 sine inner angles at each vertex
     c. Check expected mismatches probability:
        `E[mismatches] = (N_false⁴ · θ_tol⁵) / (2π²) · sin(θ_ij) / sin(min_angle)`
     d. If E[mismatches] > max_prob, skip (likely false match)
     e. Query pair-distance DB for ij, ik, ir distances
     f. Find unique (i,j,k,r) catalog match via hash-multimap intersection
     g. Check spectrality (cross-product dot-product sign)
     h. Validate all 6 distances within tolerance
  3. If unique pyramid found, run IdentifyRemainingStars to extend
- Complexity: O(PyramidAttempts × DB_queries)
- Cutoff: Default 1000 pyramid attempts
- Note: Starts searching near center of image for better coverage

**Algorithm: IdentifyRemainingStars**
- Source: `star-id.cpp:461-569`
- Purpose: After initial identification, identify additional stars using pair geometry
- Method:
  1. For each unidentified centroid, track the two closest identified stars
  2. Use the angle between unidentified→identified1 and unidentified→identified2 to find
     the "triangular angle" closest to 90°
  3. When a centroid has a near-90° angle to two identified stars (threshold = π/4),
     query the pair-distance DB for a unique third catalog star
  4. Check spectrality to ensure correct orientation
  5. If exactly 1 candidate found, identify and repeat
- Key insight: 90° angles give the most discriminating power for third-star identification

### 2.3 Attitude Estimation Algorithms

**Algorithm: Davenport Q-method (DQM)**
- Source: `attitude-estimators.cpp:12-120`
- Purpose: Optimal attitude from all identified stars (Wahba's problem)
- Method:
  1. Build attitude profile matrix B = Σ w_i · r_i · b_iᵀ
  2. Compute S = B + Bᵀ, σ = tr(B), Z = [B₂₃-B₃₂, B₃₁-B₁₃, B₁₂-B₂₁]ᵀ
  3. Build 4×4 K-matrix:
     K = [[σ, Zᵀ],
          [Z, S - σI]]
  4. Find the largest eigenvalue λ_max and corresponding eigenvector q
  5. q is the optimal quaternion
- Complexity: O(N + 4³) dominated by eigen-decomposition
- Libraries: Uses Eigen3 for eigen-decomposition

**Algorithm: TRIAD**
- Source: `attitude-estimators.cpp:122-157`
- Purpose: Fast 2-star attitude estimation
- Method:
  1. Select two identified stars
  2. Build orthonormal frames from star positions in both body (camera) and reference (catalog) frames
  3. Attitude = photo_frame · catalog_frameᵀ
- Frame construction: t₁ = v₁/|v₁|, t₂ = (v₁×v₂)/|v₁×v₂|, t₃ = t₁×t₂
- Complexity: O(1)
- Limitation: Only uses 2 stars, no noise averaging

**Algorithm: QUEST**
- Source: `attitude-estimators.cpp:184-247`
- Purpose: Fast optimal attitude without full eigen-decomposition
- Method:
  1. Build B matrix (same as DQM)
  2. Compute S, σ, Z (same as DQM)
  3. Compute characteristic polynomial coefficients:
     - δ = det(S), κ = tr(S⁻¹ · δ)
     - a = σ² - κ, b = σ² + |Z|²
     - c = δ + Z·S·Z, d = Z·S²·Z
  4. Find largest eigenvalue via Newton-Raphson on characteristic polynomial
  5. Solve for quaternion using closed-form equations
- Complexity: O(N) — no eigen-decomposition needed
- Newton-Raphson: Iterates until |f(λ)/f'(λ)| < 0.0001

### 2.4 Database Structures

**Algorithm: K-Vector Index**
- Source: `databases.cpp:64-164`
- Purpose: Accelerate range queries on sorted 1D data
- Structure:
  - Pre-compute bin boundaries for numBins equally-spaced intervals
  - Each bin stores cumulative count of values ≤ bin_boundary
  - Query: binary search → direct array access O(1) lookup
- Query: Given [q_min, q_max], return [bins[bin(q_min)-1], bins[bin(q_max)]]
- Memory: bins array of length (numBins + 1)

**Algorithm: PairDistanceKVectorDatabase**
- Source: `databases.cpp:17-295`
- Purpose: Store all inter-star distances for fast lookup
- Construction:
  1. Compute all O(N²) inter-star distances (catalog size N)
  2. Filter to [min_distance, max_distance]
  3. Sort by distance
  4. Build K-vector index over distances
  5. Store pairs as flat int16_t array: [a₁, b₁, a₂, b₂, ...]
- Query: K-vector + optional cosine refinement for exact match

## 3. Dependency Graph

```
                    ┌─────────────┐
                    │  decimal.hpp │  (typedef double|float)
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │ math_utils  │  (Vec2, Vec3, Mat3, Quat, Attitude)
                    └──┬───┬───┬──┘
                       │   │   │
          ┌────────────┘   │   └──────────┐
          ▼                ▼              ▼
    ┌──────────┐   ┌───────────┐   ┌──────────┐
    │ camera   │   │ star-utils│   │database  │
    └────┬─────┘   └─────┬─────┘   └────┬─────┘
         │               │              │
         ▼               ▼              ▼
    ┌──────────┐   ┌───────────┐   ┌──────────┐
    │centroid  │   │  star-id  │   │  KVector  │
    │ ing      │   │           │   │  Index    │
    └────┬─────┘   └─────┬─────┘   └──────────┘
         │               │
         ▼               ▼
    ┌──────────┐   ┌───────────┐
    │ pipeline │◄──│ attitude  │
    │          │   │ estimators│
    └──────────┘   └───────────┘
```

## 4. Data Flow Graph

```
INPUT IMAGE (PNG file or generated)
    │
    ▼
┌─────────────────────────────────────────────┐
│  STAGE 1: Centroid Detection                │
│  Algorithm: CoG / IWCoG / Dummy             │
│  Input: (uint8 array, width, height)        │
│  Output: list[Star] (position, radiusX,     │
│           radiusY, magnitude)               │
│  Config: centroid_algo, centroid_mag_filter,│
│          centroid_filter_brightest           │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  STAGE 2: Centroid Filtering               │
│  Algorithm: magnitude / brightest-N         │
│  Input: list[Star]                          │
│  Output: list[Star] (filtered)              │
│  Config: centroid_mag_filter (keep stars    │
│          with magnitude >= threshold),      │
│          centroid_filter_brightest (keep    │
│          top N brightest)                   │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  STAGE 3: Database Loading                 │
│  Algorithm: MultiDatabase parsing           │
│  Input: database binary file path           │
│  Output: Catalog, PairDistanceKVectorDB     │
│  Config: database_path                      │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  STAGE 4: Star Identification              │
│  Algorithm: Pyramid / Geometric Voting      │
│  Input: list[Star], Catalog, Camera,        │
│         PairDistanceKVectorDB               │
│  Output: list[StarIdentifier] (starIndex →  │
│           catalogIndex pairs)               │
│  Config: star_id_algo, angular_tolerance,   │
│          false_stars_estimate,              │
│          max_mismatch_probability           │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  STAGE 5: Attitude Estimation              │
│  Algorithm: DQM / TRIAD / QUEST             │
│  Input: list[StarIdentifier], Catalog,      │
│         list[Star], Camera                  │
│  Output: Attitude (quaternion + Euler       │
│           angles: RA, Dec, Roll)            │
│  Config: attitude_algo                      │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
            FINAL ATTITUDE SOLUTION
            - Quaternion (real, i, j, k)
            - Euler angles (RA, Dec, Roll in radians)
            - Direction Cosine Matrix (3×3)
```

## 5. Configuration System

All pipeline options (from `pipeline-options.hpp`):
```
Option                    Type      Default    Description
────────────────────────────────────────────────────────────
camera:
  png                     str       ""          PNG input path
  focal_length            float     0           Lens focal length (mm)
  pixel_size              float     -1          Pixel size (μm)
  fov                     float     20          Field of view (deg)

centroiding:
  centroid_algo           str       "cog"       "cog"/"iwcog"/"dummy"
  centroid_dummy_stars    int       5           Dummy star count
  centroid_mag_filter     float     -1          Min magnitude to keep
  centroid_filter_brightest int     -1          Keep N brightest

database:
  database_path            str       ""         Database file path

star identification:
  star_id_algo            str       "pyramid"   "py"/"gv"/"dummy"
  angular_tolerance       float     0.04        Tolerance (deg)
  false_stars_estimate    int       500         Estimated false stars
  max_mismatch_probability float    0.001       Max mismatch prob

attitude estimation:
  attitude_algo           str       "dqm"       "dqm"/"triad"/"quest"
```

Database options (from `database-options.hpp`):
```
Option                    Type      Default    Description
────────────────────────────────────────────────────────────
  min_mag                 float     100         Min catalog magnitude
  max_stars               int       10000       Max catalog stars
  min_separation          float     0.08        Min star separation (deg)
  kvector                 bool      false       Build K-vector DB
  kvector_min_distance    float     0.5         Min pair distance (deg)
  kvector_max_distance    float     15          Max pair distance (deg)
  kvector_distance_bins   int       10000       Number of K-vector bins
  output                  str       "-"         Output path
```

## 6. Porting Strategy

### Phase 1: Core Math & Types (math_utils.py)
Port Vec2, Vec3, Mat3, Quaternion, Attitude, EulerAngles, all coordinate transforms.
These are pure math with no dependencies beyond numpy.

### Phase 2: Camera Model (camera.py)
Port Camera class with spatial↔pixel transforms.
Depends on math_utils.

### Phase 3: Star Types & Catalog (star_id.py types + database.py)
Port CatalogStar, Star, StarIdentifier, Catalog narrow/filter functions.
Depends on math_utils.

### Phase 4: Database Structures (database.py)
Port KVectorIndex, PairDistanceKVectorDatabase, MultiDatabase.
Implement in-memory versions (skip binary serialization; use Python objects).
Depends on math_utils, star types.

### Phase 5: Centroiding (centroiding.py)
Port all threshold functions, CoG, IWCoG.
Depends on math_utils.

### Phase 6: Star Identification (star_id.py)
Port Geometric Voting, Pyramid, IdentifyRemainingStars.
Depends on math_utils, camera, database, star types.

### Phase 7: Attitude Estimation (attitude.py)
Port DQM, TRIAD, QUEST.
Depends on math_utils, camera, star types (for Vec3/Quaternion).

### Phase 8: Pipeline Integration (pipeline.py)
Port Pipeline class, SetPipeline factory, PipelineOutput.
Wire all stages together.

### Phase 9: Configuration (config.py)
All dataclasses for pipeline and database options.

### Phase 10: Validation (validation.py)
Comparison scripts, test matrix generation, report generation.

### Phase 11: Execution Graph (execution_graph.py)
Document each pipeline stage's IO and parameters.

## 7. Key Porting Decisions

1. **No binary serialization**: Python will use native data structures.
   - Catalog: list[CatalogStar]
   - KVectorIndex: Python object with dict/list attributes
   - PairDistanceKVectorDB: Python object with pre-computed pairs

2. **No Eigen3 dependency**: Use numpy for matrix operations where needed, but
   implement explicit math for educational transparency.

3. **Float precision**: Use Python's native float (double) throughout.

4. **No Cairo**: Skip image generation, plotting, visualization.
   Only port the analytical pipeline.

5. **No binary database format**: Read catalog from TSV or JSON.
   Build K-vector in memory at pipeline startup.

6. **Trace system**: Every stage writes JSON traces for debugging.

7. **Validation framework**: Shell scripts to run C++ LOST, capture outputs,
   compare against Python outputs.

## 8. Complete File Listing

```
lost_python/
├── README.md              - Comprehensive documentation
├── LOST_ARCHITECTURE.md   - This file
├── MIGRATION_LOG.md       - Porting progress log
├── pyproject.toml         - Python project configuration
├── config.py              - All configuration dataclasses
├── pipeline.py            - Pipeline orchestration & execution
├── execution_graph.py     - Explicit pipeline stage documentation
├── camera.py              - Camera model & coordinate transforms
├── math_utils.py          - Vector math, quaternions, attitudes
├── centroiding.py         - All centroid detection algorithms
├── star_id.py             - All star identification algorithms
├── attitude.py            - All attitude estimation algorithms
├── database.py            - Catalog loading, K-vector, DB structures
├── trace.py               - Trace/debug utilities
├── validation.py          - C++ vs Python comparison framework
├── tests/
│   ├── test_math_utils.py
│   ├── test_camera.py
│   ├── test_centroiding.py
│   ├── test_database.py
│   ├── test_star_id.py
│   ├── test_attitude.py
│   └── test_pipeline.py
├── traces/                - Runtime trace outputs
└── validation/            - Validation scripts & reports
```

## 9. Data Structure Equivalents

```
C++ Type                  Python Type
─────────────────────────────────────
decimal                   float
Vec2                      Vec2 (namedtuple or dataclass)
Vec3                      Vec3 (namedtuple or dataclass)
Mat3                      Mat3 (9-element tuple or list)
Quaternion                Quaternion (4-element namedtuple)
Star                      Star (dataclass)
CatalogStar               CatalogStar (dataclass)
StarIdentifier            StarIdentifier (dataclass)
Stars                     list[Star]
Catalog                   list[CatalogStar]
StarIdentifiers           list[StarIdentifier]
Camera                    Camera (dataclass)
Attitude                  Attitude (dataclass/quaternion)
EulerAngles               EulerAngles (dataclass)
PipelineOptions           PipelineConfig (dataclass)
DatabaseOptions           DatabaseConfig (dataclass)
PipelineOutput            PipelineOutput (dataclass)
PipelineInput             PipelineInput (protocol)
KVectorIndex              KVectorIndex (class)
PairDistanceKVectorDB     PairDistanceKVectorDB (class)
MultiDatabase             MultiDatabase (class)
```
