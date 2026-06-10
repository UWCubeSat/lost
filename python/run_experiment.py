"""
Experiment runner: Python LOST pipeline on frame_04 and frame_23.
Exact config per user spec:
  camera: sensor=1014x760, pixel_size=6.2um, focal_length=25mm
  centroid: algo=cog, mag_filter=5
  star_id: algo=pyramid, angular_tolerance=0.1, false_estimate=500, max_mismatch=0.001
  attitude: algo=dqm
"""
import math, pickle, sys, os
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from camera import Camera
from config import PipelineConfig, DatabaseConfig
from database import load_catalog, narrow_catalog, PairDistanceKVectorDatabase
from math_utils import deg_to_rad
from pipeline import Pipeline

CATALOG_PATH = "/mnt/d/WSL/py_lost/lost/bright-star-catalog.tsv"
DB_PATH = "/tmp/opencode/lost_experiment/database.pkl"

# ── Load database ──
print("Loading database...")
with open(DB_PATH, "rb") as f:
    db_data = pickle.load(f)
pair_db = db_data["pair_db"]
catalog = db_data["catalog"]
print(f"Catalog: {len(catalog)} stars, Pairs: {pair_db.num_pairs}")

def run_pipeline(image_path, label):
    print(f"\n{'='*60}")
    print(f"Python: {label}")
    print(f"{'='*60}")

    img = Image.open(image_path).convert("L")
    w, h = img.size
    print(f"Image: {w}x{h}")

    # Match C++ grayscale conversion: round(R*0.21+G*0.71+B*0.07) = round(gray*0.99)
    raw_pixels = list(img.tobytes())
    cpp_pixels = bytes(round(px * 0.99) for px in raw_pixels)
    img_bytes = cpp_pixels

    focal_pixels = 25 * 1000 / 6.2
    camera = Camera.from_center(focal_pixels, w, h)
    print(f"Focal length: {focal_pixels:.2f} px, FOV: {math.degrees(camera.fov):.4f} deg")

    cfg = PipelineConfig(
        focal_length=25,
        pixel_size=6.2,
        centroid_algo="cog",
        centroid_mag_filter=5,
        star_id_algo="pyramid",
        angular_tolerance=0.1,
        false_stars_estimate=500,
        max_mismatch_probability=0.001,
        attitude_algo="dqm",
    )

    pipe = Pipeline.from_config(cfg, catalog=catalog, pair_db=pair_db)
    output = pipe.run(img_bytes, w, h, camera=camera, trace=False)

    stars = output.stars or []
    star_ids = output.star_ids or []
    attitude = output.attitude
    num_centroids = len(stars)
    num_identified = len(star_ids)

    print(f"Centroids: {num_centroids}, Identified: {num_identified}")

    id_map = {}
    for si in star_ids:
        if si.catalog_index < len(catalog):
            id_map[si.star_index] = catalog[si.catalog_index].name

    # ── centroids.txt ──
    outdir = "/tmp/opencode/lost_experiment"
    with open(f"{outdir}/{label}_centroids.txt", "w") as f:
        f.write(f"num_actual_centroids {num_centroids}\n")
        for i, s in enumerate(stars):
            f.write(f"actual_centroid_{i}_x {s.x:.6g}\n")
            f.write(f"actual_centroid_{i}_y {s.y:.6g}\n")
            if i in id_map:
                f.write(f"actual_centroid_{i}_id {id_map[i]}\n")

    # ── attitude.txt ──
    with open(f"{outdir}/{label}_attitude.txt", "w") as f:
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

    # ── Plots ──
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 12)
    except:
        font = font_small = ImageFont.load_default()

    def draw_centroid(draw, x, y, rx, ry, color, width=2):
        r = max(rx, ry, 1)
        draw.rectangle([x-r, y-r, x+r, y+r], outline=color, width=width)

    # plot_input
    img.save(f"{outdir}/{label}_plot_input.png")

    # plot_output
    rgb = img.convert("RGB")
    draw = ImageDraw.Draw(rgb)
    meta = f"pipeline output {num_centroids} centroids   {num_identified} identified   "
    if attitude is not None and attitude.is_known():
        q = attitude.get_quaternion()
        euler = q.to_spherical()
        meta += f"RA: {math.degrees(euler.ra):.6g}  DE: {math.degrees(euler.de):.6g}  Roll: {math.degrees(euler.roll):.6g}   "
    draw.text((3, 3), meta, fill="red", font=font)
    for i, s in enumerate(stars):
        draw_centroid(draw, s.x, s.y, s.radius_x, s.radius_y, (255, 0, 0))
        if i in id_map:
            tx = s.x + max(s.radius_x, 1) + 3
            ty = s.y - max(s.radius_y, 1) + 12
            draw.text((tx, ty), str(id_map[i]), fill=(255, 100, 100), font=font_small)
    rgb.save(f"{outdir}/{label}_plot_output.png")

    # plot_centroid_indices
    rgb2 = img.convert("RGB")
    draw2 = ImageDraw.Draw(rgb2)
    draw2.text((3, 3), f"centroid indices (input) {num_centroids} centroids", fill=(255, 128, 0), font=font)
    for i, s in enumerate(stars):
        draw_centroid(draw2, s.x, s.y, s.radius_x, s.radius_y, (255, 128, 0))
        tx = s.x + max(s.radius_x, 1) + 3
        ty = s.y - max(s.radius_y, 1) + 12
        draw2.text((tx, ty), str(i), fill=(255, 128, 0), font=font_small)
    rgb2.save(f"{outdir}/{label}_plot_centroid_indices.png")

    # Print summary
    print(f"\n--- Summary for {label} ---")
    print(f"Centroids: {num_centroids}")
    print(f"Identified stars: {num_identified}")
    print(f"Catalog IDs: {[id_map.get(si.star_index, '?') for si in star_ids]}")
    unnamed = [i for i in range(num_centroids) if i not in id_map]
    if unnamed:
        print(f"Failed centroids (unidentified): {unnamed}")
    else:
        print(f"Failed centroids: none")
    if attitude is not None and attitude.is_known():
        q = attitude.get_quaternion()
        euler = q.to_spherical()
        print(f"RA: {math.degrees(euler.ra):.4f}  DE: {math.degrees(euler.de):.4f}  Roll: {math.degrees(euler.roll):.4f}")

    return stars, star_ids, attitude, id_map

# ── Run both frames ──
py04 = run_pipeline("/mnt/d/WSL/py_lost/frame_00004.png", "py_frame04")
py23 = run_pipeline("/mnt/d/WSL/py_lost/frame_00023.png", "py_frame23")
