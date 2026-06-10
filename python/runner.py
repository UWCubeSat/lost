"""
Complete runner script for the Python LOST port.
Matches C++ output format exactly for fair comparison.
"""
import math
import pickle
import sys
from PIL import Image, ImageDraw, ImageFont

from camera import Camera
from config import PipelineConfig, DatabaseConfig
from database import load_catalog, narrow_catalog, PairDistanceKVectorDatabase
from math_utils import deg_to_rad
from pipeline import Pipeline


# ── Configuration (mirrors the C++ command) ──────────────────────────

CATALOG_PATH = "../lost/bright-star-catalog.tsv"
INPUT_IMAGE = "/mnt/d/WSL/py_lost/frame_00004.png"
DB_PATH = "my-database.pkl"

pipeline_cfg = PipelineConfig(
    focal_length=25,
    pixel_size=6.2,
    centroid_algo="cog",
    centroid_mag_filter=5,
    star_id_algo="pyramid",
    angular_tolerance=0.08,
    false_stars_estimate=500,
    max_mismatch_probability=0.001,
    attitude_algo="dqm",
)

db_cfg = DatabaseConfig(
    max_stars=10000,
    min_mag=6.5,
    min_separation=0.05,
    kvector=True,
    kvector_min_distance=0.5,
    kvector_max_distance=15,
    kvector_distance_bins=10000,
)


# ── 1. Build / Load Database ─────────────────────────────────────────

print("Loading catalog...")
catalog = load_catalog(CATALOG_PATH)
narrowed = narrow_catalog(catalog, int(db_cfg.min_mag * 100), db_cfg.max_stars, deg_to_rad(db_cfg.min_separation))

print(f"Building K-vector pair database ({db_cfg.kvector_max_distance} deg max)...")
pair_db = PairDistanceKVectorDatabase(
    narrowed,
    deg_to_rad(db_cfg.kvector_min_distance),
    deg_to_rad(db_cfg.kvector_max_distance),
    db_cfg.kvector_distance_bins,
)

with open(DB_PATH, "wb") as f:
    pickle.dump({"catalog": narrowed, "pair_db": pair_db}, f)
print(f"Database saved to {DB_PATH}")


# ── 2. Load Image (matching C++ Cairo SurfaceToGrayscaleImage) ──────

print(f"Loading image: {INPUT_IMAGE}")
img = Image.open(INPUT_IMAGE).convert("L")
w, h = img.size
print(f"Image size: {w}x{h}")

# C++ SurfaceToGrayscaleImage (io.cpp:144-151) converts ARGB32 to grayscale
# using: round(R*0.21 + G*0.71 + B*0.07).
# For a grayscale PNG, Cairo sets R=G=B=gray, so the result is
# round(gray * 0.99) — losing 1% of brightness since weights sum to 0.99.
# Pillow's convert("L") preserves exact gray values.
# We replicate the C++ conversion for exact numerical match:
raw_pixels = list(img.tobytes())
cpp_pixels = bytes(round(px * 0.99) for px in raw_pixels)
img_bytes = cpp_pixels

focal_pixels = pipeline_cfg.focal_length * 1000 / pipeline_cfg.pixel_size
camera = Camera.from_center(focal_pixels, w, h)


# ── 3. Run Pipeline ──────────────────────────────────────────────────

print("Running pipeline...")
pipeline = Pipeline.from_config(pipeline_cfg, catalog=narrowed, pair_db=pair_db)
output = pipeline.run(img_bytes, w, h, camera=camera, trace=False)

stars = output.stars or []
star_ids = output.star_ids or []
attitude = output.attitude


# ── 4. Write centroids.txt (C++ format) ──────────────────────────────

print("\nWriting py_centroids.txt...")
with open("py_centroids.txt", "w") as f:
    f.write(f"num_actual_centroids {len(stars)}\n")
    for i, s in enumerate(stars):
        f.write(f"actual_centroid_{i}_x {s.x:.6g}\n")
        f.write(f"actual_centroid_{i}_y {s.y:.6g}\n")
        # Check if this centroid was identified
        for si in star_ids:
            if si.star_index == i:
                name = narrowed[si.catalog_index].name if si.catalog_index < len(narrowed) else -1
                f.write(f"actual_centroid_{i}_id {name}\n")
                break


# ── 5. Write attitude.txt (C++ format) ───────────────────────────────

print("Writing py_attitude.txt...")
with open("py_attitude.txt", "w") as f:
    if attitude is not None and attitude.is_known():
        f.write("attitude_known 1\n")
        q = attitude.get_quaternion()
        euler = q.to_spherical()
        f.write(f"attitude_ra {math.degrees(euler.ra):.6g}\n")
        f.write(f"attitude_de {math.degrees(euler.de):.6g}\n")
        f.write(f"attitude_roll {math.degrees(euler.roll):.6g}\n")
        f.write(f"attitude_i {q.i:.6g}\n")
        f.write(f"attitude_j {q.j:.6g}\n")
        f.write(f"attitude_k {q.k:.6g}\n")
        f.write(f"attitude_real {q.real:.6g}\n")
    else:
        f.write("attitude_known 0\n")


# ── 6. Generate plot images (match C++ SurfacePlot style) ────────────


def draw_centroid(draw, x, y, radius_x, radius_y, color, width=2):
    """Draw a rectangle around a centroid (matching C++ cairo_rectangle)."""
    r = max(radius_x, radius_y, 1)
    x1 = x - r
    y1 = y - r
    x2 = x + r
    y2 = y + r
    draw.rectangle([x1, y1, x2, y2], outline=color, width=width)


print("Writing py_plot_input.png...")
# --plot-raw-input: just the raw image, no annotations
img.save("py_plot_input.png")


print("Writing py_plot_output.png...")
# --plot-output: red annotations, star IDs (HR names), attitude text
rgb = img.convert("RGB")
draw = ImageDraw.Draw(rgb)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
    font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 12)
except (IOError, OSError):
    font = ImageFont.load_default()
    font_small = font

# Metadata string (matching C++ SurfacePlot)
identified_set = {si.star_index for si in star_ids}
metadata = f"pipeline output {len(stars)} centroids   {len(star_ids)} identified   "
if attitude is not None and attitude.is_known():
    q = attitude.get_quaternion()
    euler = q.to_spherical()
    metadata += f"RA: {math.degrees(euler.ra):.6g}  DE: {math.degrees(euler.de):.6g}  Roll: {math.degrees(euler.roll):.6g}   "

draw.text((3, 3), metadata, fill="red", font=font)

# Draw centroids (red, alpha 0.5 simulated via stipple)
for i, s in enumerate(stars):
    draw_centroid(draw, s.x, s.y, s.radius_x, s.radius_y, (255, 0, 0), width=2)

    # Draw identified star names next to centroid
    for si in star_ids:
        if si.star_index == i:
            name = narrowed[si.catalog_index].name if si.catalog_index < len(narrowed) else -1
            text_x = s.x + max(s.radius_x, 1) + 3
            text_y = s.y - max(s.radius_y, 1) + 12
            draw.text((text_x, text_y), str(name), fill=(255, 100, 100), font=font_small)
            break

rgb.save("py_plot_output.png")


print("Writing py_plot_centroid_indices.png...")
# C++ style: orange annotations, centroid index numbers, no attitude text
rgb2 = img.convert("RGB")
draw2 = ImageDraw.Draw(rgb2)

metadata2 = f"centroid indices (input) {len(stars)} centroids   "
draw2.text((3, 3), metadata2, fill=(255, 128, 0), font=font)

for i, s in enumerate(stars):
    draw_centroid(draw2, s.x, s.y, s.radius_x, s.radius_y, (255, 128, 0), width=2)
    text_x = s.x + max(s.radius_x, 1) + 3
    text_y = s.y - max(s.radius_y, 1) + 12
    draw2.text((text_x, text_y), str(i), fill=(255, 128, 0), font=font_small)

rgb2.save("py_plot_centroid_indices.png")


print("\nDone. Output files:")
for fname in ["py_centroids.txt", "py_attitude.txt", "py_plot_input.png",
              "py_plot_output.png", "py_plot_centroid_indices.png"]:
    print(f"  {fname}")
