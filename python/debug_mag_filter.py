"""
Debug script for centroid_mag_filter regression test.
Runs pipeline with mag_filter=6 and mag_filter=5, compares results.
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


CATALOG_PATH = "../lost/bright-star-catalog.tsv"
INPUT_IMAGE = "/mnt/d/WSL/py_lost/frame_00004.png"
DB_PATH = "my-database.pkl"

# ---- Load database ----
print("Loading catalog...")
catalog = load_catalog(CATALOG_PATH)
narrowed = narrow_catalog(catalog, 650, 10000, deg_to_rad(0.05))

print("Loading database...")
with open(DB_PATH, "rb") as f:
    db_data = pickle.load(f)
pair_db = db_data["pair_db"]
narrowed = db_data["catalog"]

# ---- Load image ----
print(f"Loading image: {INPUT_IMAGE}")
img = Image.open(INPUT_IMAGE).convert("L")
w, h = img.size
print(f"Image size: {w}x{h}")
raw_pixels = list(img.tobytes())
cpp_pixels = bytes(round(px * 0.99) for px in raw_pixels)
img_bytes = cpp_pixels

# ---- Camera (same for both cases) ----
focal_length_mm = 25
pixel_size_um = 6.2
focal_pixels = focal_length_mm * 1000 / pixel_size_um
camera = Camera.from_center(focal_pixels, w, h)


def run_case(mag_filter: int, label: str):
    """Run pipeline with given mag_filter, return results with full debug."""
    cfg = PipelineConfig(
        focal_length=focal_length_mm,
        pixel_size=pixel_size_um,
        centroid_algo="cog",
        centroid_mag_filter=mag_filter,
        star_id_algo="pyramid",
        angular_tolerance=0.08,
        false_stars_estimate=500,
        max_mismatch_probability=0.001,
        attitude_algo="dqm",
    )
    pipe = Pipeline.from_config(cfg, catalog=narrowed, pair_db=pair_db)
    output = pipe.run(img_bytes, w, h, camera=camera, trace=False)

    stars = output.stars or []
    star_ids = output.star_ids or []
    attitude = output.attitude

    # Build lookup: star_index -> identified info
    identified_map = {}
    for si in star_ids:
        identified_map[si.star_index] = {
            "catalog_index": si.catalog_index,
            "name": narrowed[si.catalog_index].name if si.catalog_index < len(narrowed) else -1,
        }

    return {
        "label": label,
        "stars": stars,
        "star_ids": star_ids,
        "attitude": attitude,
        "identified_map": identified_map,
        "num_centroids": len(stars),
        "num_identified": len(star_ids),
    }


# ---- Run both cases ----
case6 = run_case(6, "mag_filter=6")
case5 = run_case(5, "mag_filter=5")

# ---- Comparison ----
print("\n" + "=" * 80)
print("COMPARISON: CASE A (mag=6) vs CASE B (mag=5)")
print("=" * 80)

print(f"\nCASE A: {case6['num_centroids']} centroids, {case6['num_identified']} identified")
print(f"CASE B: {case5['num_centroids']} centroids, {case5['num_identified']} identified")

# Print full centroid table
print("\n--- Centroid Comparison ---")
print(f"{'idx':>4} {'x_A':>10} {'y_A':>10} {'x_B':>10} {'y_B':>10} {'mag_A':>6} {'mag_B':>6} {'ID_A':>8} {'ID_B':>8} {'match':>6}")
print("-" * 80)

# Match centroids by proximity (same physical star)
matched_B = set()
for i_A, s_A in enumerate(case6["stars"]):
    id_A = case6["identified_map"].get(i_A, {}).get("name", "---")
    best_dist = 999
    best_j = -1
    for j_B, s_B in enumerate(case5["stars"]):
        dist = math.hypot(s_A.x - s_B.x, s_A.y - s_B.y)
        if dist < best_dist:
            best_dist = dist
            best_j = j_B
    if best_j >= 0 and best_dist < 2.0:
        id_B = case5["identified_map"].get(best_j, {}).get("name", "---")
        matched_B.add(best_j)
        match_str = "OK" if best_dist < 1 else "NEAR"
        print(f"{i_A:>4} {s_A.x:>10.4f} {s_A.y:>10.4f} {case5['stars'][best_j].x:>10.4f} {case5['stars'][best_j].y:>10.4f} {s_A.magnitude:>6} {case5['stars'][best_j].magnitude:>6} {str(id_A):>8} {str(id_B):>8} {match_str:>6}")
    else:
        print(f"{i_A:>4} {s_A.x:>10.4f} {s_A.y:>10.4f} {'---':>10} {'---':>10} {s_A.magnitude:>6} {'---':>6} {str(id_A):>8} {'---':>8} {'MISS':>6}")

# Print unmatched B centroids
for j_B, s_B in enumerate(case5["stars"]):
    if j_B not in matched_B:
        id_B = case5["identified_map"].get(j_B, {}).get("name", "---")
        print(f"{'NEW':>4} {'---':>10} {'---':>10} {s_B.x:>10.4f} {s_B.y:>10.4f} {'---':>6} {s_B.magnitude:>6} {'---':>8} {str(id_B):>8} {'NEW':>6}")

# ---- Stability analysis ----
print("\n--- Identification Stability ---")
print(f"{'idx_A':>6} {'idx_B':>6} {'x':>10} {'y':>10} {'mag_A':>6} {'mag_B':>6} {'id_A':>8} {'id_B':>8} {'status':>12}")
print("-" * 72)

unstable = []
for i_A, s_A in enumerate(case6["stars"]):
    id_A = case6["identified_map"].get(i_A, {}).get("name", "---")
    best_dist = 999
    best_j = -1
    for j_B, s_B in enumerate(case5["stars"]):
        dist = math.hypot(s_A.x - s_B.x, s_A.y - s_B.y)
        if dist < best_dist:
            best_dist = dist
            best_j = j_B
    if best_j >= 0 and best_dist < 2.0:
        id_B = case5["identified_map"].get(best_j, {}).get("name", "---")
        status = "STABLE-ID" if (id_A == id_B and id_A != "---") else "UNSTABLE" if (id_A != "---" and (id_B == "---" or id_A != id_B)) else "BOTH-UNK"
        if id_A != "---" and id_B == "---":
            unstable.append((i_A, best_j, s_A.x, s_A.y, id_A, id_B))
        print(f"{i_A:>6} {best_j:>6} {s_A.x:>10.4f} {s_A.y:>10.4f} {s_A.magnitude:>6} {case5['stars'][best_j].magnitude:>6} {str(id_A):>8} {str(id_B):>8} {status:>12}")
    else:
        id_B_val = f"idx {[j for j in range(len(case5['stars'])) if math.hypot(s_A.x - case5['stars'][j].x, s_A.y - case5['stars'][j].y) < 3]}"
        print(f"{i_A:>6} {'---':>6} {s_A.x:>10.4f} {s_A.y:>10.4f} {s_A.magnitude:>6} {'---':>6} {str(id_A):>8} {'---':>8} {'MISS':>12}")

if unstable:
    print(f"\n!!! {len(unstable)} stars IDENTIFIED in CASE A but UNIDENTIFIED in CASE B:")
    for (i_A, j_B, x, y, id_A, id_B) in unstable:
        print(f"  Centroid A[{i_A}] (B[{j_B}]) at ({x:.2f}, {y:.2f}): Star {id_A} lost!")

# ---- Root cause analysis ----
print("\n--- Root Cause Analysis ---")
print("Checking pyramid search star selection...")
stars6 = case6["stars"]
stars5 = case5["stars"]
print(f"  mag=6: first 10 star indices have magnitudes:")
for i in range(min(10, len(stars6))):
    s = stars6[i]
    id_info = case6["identified_map"].get(i, {})
    name = id_info.get("name", "---")
    print(f"    idx={i}: x={s.x:.2f}, y={s.y:.2f}, mag={s.magnitude}, id={name}")

print(f"\n  mag=5: first 10 star indices have magnitudes:")
for i in range(min(10, len(stars5))):
    s = stars5[i]
    id_info = case5["identified_map"].get(i, {})
    name = id_info.get("name", "---")
    print(f"    idx={i}: x={s.x:.2f}, y={s.y:.2f}, mag={s.magnitude}, id={name}")

# Show which stars appear in first PYRAMID_LIMIT for each case
print(f"\n  PYRAMID_LIMIT = 10")
print(f"  mag=6 stars in pyramid range (indices 0-9):")
for i in range(min(10, len(stars6))):
    s = stars6[i]
    print(f"    [{i}] mag={s.magnitude} @ ({s.x:.2f}, {s.y:.2f})")

print(f"\n  mag=5 stars in pyramid range (indices 0-9):")
for i in range(min(10, len(stars5))):
    s = stars5[i]
    print(f"    [{i}] mag={s.magnitude} @ ({s.x:.2f}, {s.y:.2f})")

# Check if additional faint stars in mag=5 push some brighter stars out of pyramid range
brightest6 = sorted(stars6, key=lambda s: s.magnitude, reverse=True)
brightest5 = sorted(stars5, key=lambda s: s.magnitude, reverse=True)
print(f"\n  Brightest 10 stars (mag=6):")
for i in range(min(10, len(brightest6))):
    s = brightest6[i]
    print(f"    mag={s.magnitude} @ ({s.x:.2f}, {s.y:.2f})")

print(f"\n  Brightest 10 stars (mag=5):")
for i in range(min(10, len(brightest5))):
    s = brightest5[i]
    print(f"    mag={s.magnitude} @ ({s.x:.2f}, {s.y:.2f})")

# Check if the pyramid search star set differs
stars6_set = set((round(s.x, 2), round(s.y, 2)) for s in stars6[:10])
stars5_set = set((round(s.x, 2), round(s.y, 2)) for s in stars5[:10])
print(f"\n  First-10 star coordinates (rounded to 2dp):")
print(f"    mag=6: {sorted(stars6_set)}")
print(f"    mag=5: {sorted(stars5_set)}")
print(f"    Same? {'YES' if stars6_set == stars5_set else 'NO'}")
if stars6_set != stars5_set:
    print(f"    Only in mag=6: {stars6_set - stars5_set}")
    print(f"    Only in mag=5: {stars5_set - stars6_set}")

# ---- Output centroids, attitude, plots ----
def draw_centroid(draw, x, y, radius_x, radius_y, color, width=2):
    r = max(radius_x, radius_y, 1)
    x1 = x - r
    y1 = y - r
    x2 = x + r
    y2 = y + r
    draw.rectangle([x1, y1, x2, y2], outline=color, width=width)

for case, suffix in [(case6, "mag6"), (case5, "mag5")]:
    stars = case["stars"]
    star_ids = case["star_ids"]
    attitude = case["attitude"]
    identified_map = case["identified_map"]

    # centroids.txt
    with open(f"centroids_{suffix}.txt", "w") as f:
        f.write(f"num_actual_centroids {len(stars)}\n")
        for i, s in enumerate(stars):
            f.write(f"actual_centroid_{i}_x {s.x:.6g}\n")
            f.write(f"actual_centroid_{i}_y {s.y:.6g}\n")
            if i in identified_map:
                name = identified_map[i]["name"]
                f.write(f"actual_centroid_{i}_id {name}\n")

    # attitude.txt
    with open(f"attitude_{suffix}.txt", "w") as f:
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

    # plot images
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 12)
    except:
        font = ImageFont.load_default()
        font_small = font

    # plot_input - raw image
    img.save(f"plot_input_{suffix}.png")

    # plot_output
    rgb = img.convert("RGB")
    draw = ImageDraw.Draw(rgb)
    identified_set = {si.star_index for si in star_ids}
    metadata = f"pipeline output {len(stars)} centroids   {len(star_ids)} identified   "
    if attitude is not None and attitude.is_known():
        q = attitude.get_quaternion()
        euler = q.to_spherical()
        metadata += f"RA: {math.degrees(euler.ra):.6g}  DE: {math.degrees(euler.de):.6g}  Roll: {math.degrees(euler.roll):.6g}   "
    draw.text((3, 3), metadata, fill="red", font=font)
    for i, s in enumerate(stars):
        draw_centroid(draw, s.x, s.y, s.radius_x, s.radius_y, (255, 0, 0), width=2)
        if i in identified_map:
            name = identified_map[i]["name"]
            text_x = s.x + max(s.radius_x, 1) + 3
            text_y = s.y - max(s.radius_y, 1) + 12
            draw.text((text_x, text_y), str(name), fill=(255, 100, 100), font=font_small)
    rgb.save(f"plot_output_{suffix}.png")

    # plot_centroid_indices
    rgb2 = img.convert("RGB")
    draw2 = ImageDraw.Draw(rgb2)
    metadata2 = f"centroid indices (input) {len(stars)} centroids   "
    draw2.text((3, 3), metadata2, fill=(255, 128, 0), font=font)
    for i, s in enumerate(stars):
        draw_centroid(draw2, s.x, s.y, s.radius_x, s.radius_y, (255, 128, 0), width=2)
        text_x = s.x + max(s.radius_x, 1) + 3
        text_y = s.y - max(s.radius_y, 1) + 12
        draw2.text((text_x, text_y), str(i), fill=(255, 128, 0), font=font_small)
    rgb2.save(f"plot_centroid_indices_{suffix}.png")

print("\nDone! All outputs written.")
