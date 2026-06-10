"""
====================================================
EXECUTION GRAPH — Pipeline Stage Documentation
====================================================

Purpose:
    Documents every stage of the pipeline with:
        - Input schema
        - Output schema
        - Algorithm used
        - Configuration parameters
        - Mathematical operations

    This is the definitive reference for understanding what
    happens at each step of the star-tracking process.

Stage Index:
    Stage 1: Input Image → raw pixel data
    Stage 2: Centroid Detection → (x, y) coordinates
    Stage 3: Centroid Filtering → filtered centroids
    Stage 4: Database Loading → catalog + pair distances
    Stage 5: Star Identification → star → catalog mapping
    Stage 6: Attitude Estimation → orientation quaternion
    Stage 7: Final Solution → Euler angles + quaternion
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class Port:
    """A data port (input or output) of a pipeline stage."""

    name: str
    type: str
    description: str
    shape: Optional[str] = None  # e.g., "(N,)", "(3,3)", etc.


@dataclass
class ConfigParam:
    """A configuration parameter for a pipeline stage."""

    name: str
    type: str
    default: Any
    description: str
    choices: Optional[List[str]] = None


@dataclass
class PipelineStage:
    """Description of a single pipeline stage."""

    name: str
    number: int
    purpose: str
    algorithm: str
    algorithm_variants: List[str]
    inputs: List[Port]
    outputs: List[Port]
    config: List[ConfigParam]
    mathematical_operations: List[str]
    complexity: str
    failure_modes: List[str]
    notes: str = ""


# ─────────────────────────────────────────
# STAGE DEFINITIONS
# ─────────────────────────────────────────


STAGES = [
    PipelineStage(
        name="Input Image",
        number=1,
        purpose="Read raw image data from PNG file or generated image.",
        algorithm="PNG decode (via pillow) or image generation",
        algorithm_variants=["png_file", "generated_image"],
        inputs=[
            Port("file_path", "str", "Path to PNG image file"),
            Port("camera_config", "CameraConfig", "Focal length, pixel size, FOV"),
        ],
        outputs=[
            Port("image", "bytes", "Raw grayscale pixel data, 1 byte per pixel"),
            Port("width", "int", "Image width in pixels"),
            Port("height", "int", "Image height in pixels"),
            Port("camera", "Camera", "Camera model with transforms"),
        ],
        config=[
            ConfigParam("png", "str", "", "Path to PNG input image"),
            ConfigParam("focal_length", "float", 0.0, "Camera focal length in mm"),
            ConfigParam("pixel_size", "float", -1.0, "Pixel size in microns"),
            ConfigParam("fov", "float", 20.0, "Field of view in degrees"),
        ],
        mathematical_operations=[
            "Grayscale conversion (luminosity method: 0.21R + 0.71G + 0.07B)",
            "Focal length calculation: f_px = f_mm * 1000 / pixel_size_um",
            "Focal length from FOV: f = W / (2 * tan(FOV/2))",
        ],
        complexity="O(W × H)",
        failure_modes=[
            "File not found or corrupt",
            "Incorrect camera parameters lead to failed star ID downstream",
        ],
    ),
    PipelineStage(
        name="Centroid Detection",
        number=2,
        purpose="Detect bright star centroids in the grayscale image.",
        algorithm="Center of Gravity (CoG) — flood-fill + intensity-weighted average",
        algorithm_variants=["cog", "iwcog", "dummy"],
        inputs=[
            Port("image", "bytes", "Grayscale pixel data, 1 byte per pixel"),
            Port("width", "int", "Image width"),
            Port("height", "int", "Image height"),
        ],
        outputs=[
            Port("centroids", "List[Star]", "Detected star centroids with pixel coordinates"),
        ],
        config=[
            ConfigParam("centroid_algo", "str", "cog", "Centroid algorithm to use",
                        choices=["cog", "iwcog", "dummy"]),
            ConfigParam("centroid_dummy_stars", "int", 5,
                        "Number of dummy centroids (dummy mode only)"),
            ConfigParam("centroid_mag_filter", "float", -1.0,
                        "Minimum magnitude threshold to keep centroid"),
            ConfigParam("centroid_filter_brightest", "int", -1,
                        "Keep only N brightest centroids"),
        ],
        mathematical_operations=[
            "Threshold: T = mean(pixels) + 5 × stddev(pixels)",
            "Flood-fill: 4-connected region growing from bright seed pixels",
            "Centroid: x_c = Σ(x_i × I_i) / Σ(I_i), y_c = Σ(y_i × I_i) / Σ(I_i)",
            "IWCoG: Iterative Gaussian-weighted refinement of centroid estimate",
            "Radius: half of bounding box dimensions",
        ],
        complexity="O(W × H + N × pixels_per_star)",
        failure_modes=[
            "Too many false positives from noise",
            "Star bleeding into adjacent stars (blending)",
            "Stars too dim relative to threshold",
            "Stars on image edge are rejected",
        ],
        notes="The +0.5 pixel offset converts from corner-based to center-based coordinates.",
    ),
    PipelineStage(
        name="Centroid Filtering",
        number=3,
        purpose="Remove low-quality centroids before star identification.",
        algorithm="Magnitude threshold and/or brightest-N filter",
        algorithm_variants=["magnitude_filter", "brightest_n_filter"],
        inputs=[
            Port("centroids", "List[Star]", "Unfiltered centroids from centroiding stage"),
        ],
        outputs=[
            Port("filtered_centroids", "List[Star]", "Filtered centroids"),
        ],
        config=[
            ConfigParam("centroid_mag_filter", "float", -1.0,
                        "Keep centroids with magnitude >= this value"),
            ConfigParam("centroid_filter_brightest", "int", -1,
                        "Keep only this many brightest centroids"),
        ],
        mathematical_operations=[
            "Sort by magnitude descending",
            "Keep stars with magnitude >= max(min_mag, Nth_brightest_mag)",
        ],
        complexity="O(N log N)",
        failure_modes=[
            "Filter too aggressive → too few stars for star ID",
            "Filter too lenient → many false stars confuse star ID",
        ],
    ),
    PipelineStage(
        name="Database Loading",
        number=4,
        purpose="Load star catalog and pair-distance database for star identification.",
        algorithm="Bright Star Catalog TSV parsing + K-vector index construction",
        algorithm_variants=["binary_database", "in_memory_generation"],
        inputs=[
            Port("database_path", "str", "Path to .dat database file"),
        ],
        outputs=[
            Port("catalog", "Catalog", "List of catalog stars with positions and magnitudes"),
            Port("pair_db", "PairDistanceKVectorDatabase",
                 "Inter-star distance database with K-vector index"),
        ],
        config=[
            ConfigParam("database_path", "str", "",
                        "Path to serialized database file"),
            ConfigParam("min_mag", "float", 100.0,
                        "Maximum magnitude (×100) for catalog stars"),
            ConfigParam("max_stars", "int", 10000,
                        "Maximum number of catalog stars"),
            ConfigParam("min_separation", "float", 0.08,
                        "Minimum angular separation between catalog stars (deg)"),
        ],
        mathematical_operations=[
            "Spherical to Cartesian: (RA, Dec) → (cos(RA)cos(Dec), sin(RA)cos(Dec), sin(Dec))",
            "Inter-star distance: θ = acos(v1 · v2) for unit vectors",
            "K-vector: bin boundaries at equally-spaced intervals, cumulative counts",
            "Pair-distance DB query: O(1) range lookup via K-vector index",
        ],
        complexity="Catalog filter: O(N²) for min_separation check",
        failure_modes=[
            "Database file not found or corrupt",
            "Catalog too small for reliable star ID",
            "K-vector bins too coarse → many false positive pairs",
        ],
    ),
    PipelineStage(
        name="Star Identification",
        number=5,
        purpose="Match detected centroids to catalog stars.",
        algorithm="Pyramid — 4-star pattern matching with probabilistic validation",
        algorithm_variants=["pyramid", "geometric_voting", "dummy"],
        inputs=[
            Port("centroids", "List[Star]", "Filtered centroids from previous stage"),
            Port("catalog", "Catalog", "Star catalog with positions"),
            Port("camera", "Camera", "Camera model for pixel→spatial conversion"),
            Port("pair_db", "PairDistanceKVectorDatabase",
                 "Inter-star distance database"),
        ],
        outputs=[
            Port("star_ids", "List[StarIdentifier]",
                 "Mapping from centroid index to catalog index"),
            Port("identified_count", "int", "Number of successfully identified stars"),
        ],
        config=[
            ConfigParam("star_id_algo", "str", "pyramid",
                        "Star identification algorithm",
                        choices=["pyramid", "gv", "dummy"]),
            ConfigParam("angular_tolerance", "float", 0.04,
                        "Angular tolerance for distance matching (deg)"),
            ConfigParam("false_stars_estimate", "int", 500,
                        "Estimated number of false stars on celestial sphere"),
            ConfigParam("max_mismatch_probability", "float", 0.001,
                        "Maximum acceptable probability of mismatch"),
        ],
        mathematical_operations=[
            "Pixel→spatial: (x, y) → (1, -(x-cx)/f, -(y-cy)/f) → normalize",
            "Great-circle angle: θ = acos(v1 · v2)",
            "Expected mismatches: E[m] = N_f⁴ × θ_tol⁵ × sin(θ_ij) / (2π² × sin(min_θ_inner))",
            "Spectrality check: sign of (v1 × v2 · v3) to ensure correct handedness",
            "IdentifyRemainingStars: triangular angle near 90° → unique third star",
        ],
        complexity="Pyramid: O(N_iter × DB_queries), GV: O(N_centroids² × pairs_per_query)",
        failure_modes=[
            "Too few centroids (< 4 for pyramid, < 2 for GV)",
            "Angular tolerance too tight or too loose",
            "High false-star rate causes non-unique matches",
            "Pyramid cutoff reached without finding a match",
            "Centroid error exceeds angular tolerance",
        ],
        notes="The pyramid algorithm starts searching near the center of the image for better coverage.",
    ),
    PipelineStage(
        name="Attitude Estimation",
        number=6,
        purpose="Compute spacecraft orientation from identified stars.",
        algorithm="Davenport Q-method — K-matrix eigen-decomposition",
        algorithm_variants=["dqm", "triad", "quest"],
        inputs=[
            Port("star_ids", "List[StarIdentifier]",
                 "Identified star mappings"),
            Port("catalog", "Catalog", "Star catalog with reference positions"),
            Port("camera", "Camera", "Camera model"),
            Port("centroids", "List[Star]", "Detected centroid positions"),
        ],
        outputs=[
            Port("attitude", "Attitude", "Estimated spacecraft orientation"),
            Port("quaternion", "Quaternion",
                 "Rotation quaternion (real, i, j, k)"),
            Port("euler_angles", "EulerAngles",
                 "Right ascension, declination, roll (radians)"),
            Port("dcm", "Mat3",
                 "Direction cosine matrix (3×3 rotation matrix)"),
        ],
        config=[
            ConfigParam("attitude_algo", "str", "dqm",
                        "Attitude estimation algorithm",
                        choices=["dqm", "triad", "quest"]),
        ],
        mathematical_operations=[
            "Attitude profile matrix: B = Σ w_i × r_i × b_i^T",
            "DQM: Build K-matrix, find dominant eigenvector → optimal quaternion",
            "TRIAD: Build orthonormal frames from 2 stars, A = frame_body × frame_ref^T",
            "QUEST: Newton-Raphson eigenvalue estimation, then analytical quaternion",
            "Wahba loss: J(A) = Σ w_i × |b_i - A×r_i|²",
        ],
        complexity="DQM: O(N + 4³), TRIAD: O(1), QUEST: O(N + iterations)",
        failure_modes=[
            "Fewer than 2 identified stars (attitude unknown)",
            "Poor star geometry (all stars collinear)",
            "Incorrect star identifications cause large errors",
        ],
        notes="DQM is recommended for accuracy. TRIAD is fastest but uses only 2 stars.",
    ),
    PipelineStage(
        name="Final Solution",
        number=7,
        purpose="Produce the final attitude solution with diagnostics.",
        algorithm="Combination of all previous stages",
        algorithm_variants=[],
        inputs=[
            Port("attitude", "Attitude", "Estimated attitude"),
            Port("star_ids", "List[StarIdentifier]", "Star identifications"),
            Port("centroids", "List[Star]", "Detected centroids"),
        ],
        outputs=[
            Port("quaternion", "Quaternion",
                 "Final rotation quaternion (real, i, j, k)"),
            Port("ra", "float", "Right ascension in degrees"),
            Port("de", "float", "Declination in degrees"),
            Port("roll", "float", "Roll in degrees"),
            Port("num_identified", "int", "Number of identified stars"),
            Port("num_centroids", "int", "Number of detected centroids"),
            Port("execution_times", "Dict[str, float]",
                 "Per-stage execution times in seconds"),
        ],
        config=[],
        mathematical_operations=[
            "Quaternion → Euler angles (z-y'-x'' convention)",
            "RA = atan2(2(-qw·qz + qx·qy), 1 - 2(qy² + qz²))",
            "Dec = -asin(2(-qw·qy - qx·qz))",
            "Roll = -atan2(2(-qw·qx + qy·qz), 1 - 2(qx² + qy²))",
        ],
        complexity="O(1)",
        failure_modes=[
            "Attitude unknown (not enough stars identified)",
        ],
    ),
]


def print_stage_summary() -> None:
    """Print a summary of all pipeline stages."""
    for stage in STAGES:
        print(f"\n{'=' * 60}")
        print(f"Stage {stage.number}: {stage.name}")
        print(f"{'=' * 60}")
        print(f"  Purpose: {stage.purpose}")
        print(f"  Algorithm: {stage.algorithm}")
        print(f"  Complexity: {stage.complexity}")
        print(f"  Inputs:")
        for p in stage.inputs:
            print(f"    - {p.name} ({p.type}): {p.description}")
        print(f"  Outputs:")
        for p in stage.outputs:
            print(f"    - {p.name} ({p.type}): {p.description}")
        print(f"  Config:")
        for c in stage.config:
            print(f"    - {c.name} ({c.type}, default={c.default}): {c.description}")
