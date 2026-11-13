#!/usr/bin/env python3
"""
Tetra benchmark script.

Measures memory (db size, pattern count), latency (centroiding, star-id, attitude),
and accuracy (star-id correctness, attitude error) across many generated images.

Usage:
    python3 benchmark.py [--num-images 100] [--fov 20] [--max-angle 25] [--label baseline]
                         [--parallel 4]
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed


LOST_BIN = "./lost"
BENCH_DIR = "benchmarks"
DB_PATH_TEMPLATE = "/tmp/tetra-bench-{label}.dat"


def run(cmd, timeout=600):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=timeout
    )
    return result.stdout, result.stderr, result.returncode


def generate_database(algo, max_angle, db_path, db_flags=""):
    if algo == "tetra":
        cmd = (
            f"{LOST_BIN} database"
            f" --tetra"
            f" --tetra-max-fov {max_angle}"
            f" {db_flags}"
            f" --output {db_path}"
        )
    else:
        cmd = (
            f"{LOST_BIN} database"
            f" --kvector"
            f" --output {db_path}"
        )

    t0 = time.monotonic()
    stdout, stderr, rc = run(cmd)
    elapsed = time.monotonic() - t0

    if rc != 0:
        print(f"ERROR: database generation failed: {stderr}", file=sys.stderr)
        sys.exit(1)

    combined = stdout + stderr
    num_patterns = 0
    m = re.search(r"found (\d+) patterns", combined)
    if m:
        num_patterns = int(m.group(1))

    db_size = os.path.getsize(db_path)

    return {
        "db_generation_time_s": round(elapsed, 3),
        "db_size_bytes": db_size,
        "db_size_mb": round(db_size / 1e6, 2),
        "num_patterns": num_patterns,
    }


def pipeline_cmd(db_path, fov, num_images, seed=None, star_id_algo="tetra", extra_flags=""):
    """LOST pipeline command for N generated images."""
    seed_flag = f" --generate-seed {seed}" if seed is not None else ""
    # Only include default noise params if not overridden by extra_flags
    defaults = ""
    if "--generate-read-noise-stddev" not in extra_flags:
        defaults += " --generate-read-noise-stddev 0.05"
    if "--generate-spread-stddev" not in extra_flags:
        defaults += " --generate-spread-stddev 1"
    return (
        f"{LOST_BIN} pipeline"
        f" --generate {num_images}"
        f" --generate-random-attitudes"
        f" --generate-x-resolution 1024"
        f" --generate-y-resolution 1024"
        f" --fov {fov}"
        f" --generate-exposure 1.0"
        f"{defaults}"
        f"{seed_flag}"
        f" --database {db_path}"
        f" --centroid-algo cog"
        f" --centroid-mag-filter 3"
        f" --star-id-algo {star_id_algo}"
        f" {extra_flags}"
        f" --compare-star-ids"
        f" --compare-attitudes"
        f" --attitude-algo dqm"
        f" --print-speed"
    )


def run_batch(db_path, fov, num_images, seed=None, star_id_algo="tetra", extra_flags=""):
    """Run N generated images in a single LOST invocation."""
    stdout, stderr, rc = run(pipeline_cmd(db_path, fov, num_images, seed, star_id_algo, extra_flags))
    return stdout + stderr


def parse_results(combined, num_images, seed=None):
    """Parse all metrics from LOST output."""
    result = {"num_images": num_images}

    # Star-ID accuracy
    m = re.search(r"starid_num_images_correct (\d+)", combined)
    result["images_correct"] = int(m.group(1)) if m else 0

    m = re.search(r"starid_num_images_incorrect (\d+)", combined)
    result["images_incorrect"] = int(m.group(1)) if m else 0

    # Parse ground truth attitudes for all images
    gt_attitudes = []
    for m in re.finditer(r"(?:Random|Given) attitude:\n([\d.eE+-]+), ([\d.eE+-]+), ([\d.eE+-]+)", combined):
        gt_attitudes.append({"ra": float(m.group(1)), "de": float(m.group(2)), "roll": float(m.group(3))})

    # Parse per-image star-ID stats
    starid_by_image = {}
    for m in re.finditer(r"starid_image (\d+) correct=(\d+) incorrect=(\d+) total=(\d+)", combined):
        starid_by_image[int(m.group(1))] = {
            "correct": int(m.group(2)),
            "incorrect": int(m.group(3)),
            "total": int(m.group(4)),
        }

    # Identifiable stars per image = centroids (after mag filter) that match a real catalog star,
    # as counted by the star-ID comparator. Fewer than 4 means no 4-star pattern can be formed.
    centroid_counts = [starid_by_image[i]["total"] for i in sorted(starid_by_image)]
    result["avg_centroids"] = round(sum(centroid_counts) / len(centroid_counts), 1) if centroid_counts else 0
    result["min_centroids"] = min(centroid_counts) if centroid_counts else 0
    result["images_too_few_centroids"] = sum(1 for c in centroid_counts if c < 4)

    # Per-image attitude errors with details
    error_pattern = re.compile(
        r"attitude_error_image ([\d.eE+-]+)"
        r" image=(\d+)"
        r" expected_ra=([\d.eE+-]+)"
        r" expected_de=([\d.eE+-]+)"
        r" expected_roll=([\d.eE+-]+)"
        r" actual_ra=([\d.eE+-]+)"
        r" actual_de=([\d.eE+-]+)"
        r" actual_roll=([\d.eE+-]+)"
    )
    attitude_errors = []
    error_details = []
    images_with_attitude = set()
    for m in error_pattern.finditer(combined):
        err = float(m.group(1))
        img_idx = int(m.group(2))
        attitude_errors.append(err)
        images_with_attitude.add(img_idx)
        error_details.append({
            "error_deg": round(err, 6),
            "image": img_idx,
            "expected": {"ra": float(m.group(3)), "de": float(m.group(4)), "roll": float(m.group(5))},
            "actual": {"ra": float(m.group(6)), "de": float(m.group(7)), "roll": float(m.group(8))},
        })
    result["attitude_errors_deg"] = attitude_errors

    # Build per-image results for every image
    error_by_image = {d["image"]: d for d in error_details}
    per_image = []
    for i in range(len(gt_attitudes)):
        sid = starid_by_image.get(i, {})
        if i in error_by_image:
            d = error_by_image[i]
            err = d["error_deg"]
            has_incorrect_ids = sid.get("incorrect", 0) > 0
            if err > 0.1:
                reason = "wrong_attitude"
            elif has_incorrect_ids:
                reason = "wrong_ids_ok_attitude"
            else:
                reason = "correct"
            per_image.append({**d, "reason": reason, "starid": sid})
        else:
            reason = "too_few_centroids" if sid.get("total", 0) < 4 else "failed_to_match"
            per_image.append({
                "error_deg": None,
                "image": i,
                "expected": gt_attitudes[i],
                "actual": None,
                "reason": reason,
                "starid": sid,
            })
    for p in per_image:
        p["seed"] = seed
        p["batch_image"] = p["image"]
    result["per_image"] = per_image
    result["failures"] = [p for p in per_image if p["reason"] != "correct"]

    if attitude_errors:
        sorted_errs = sorted(attitude_errors)
        result["attitude_error_mean_deg"] = round(sum(attitude_errors) / len(attitude_errors), 6)
        result["attitude_error_median_deg"] = round(sorted_errs[len(sorted_errs) // 2], 6)
        result["attitude_error_max_deg"] = round(max(attitude_errors), 6)
        result["attitude_error_p95_deg"] = round(sorted_errs[int(0.95 * len(sorted_errs))], 6)

    # Timing from --print-speed (nanoseconds)
    for prefix in ["centroiding", "starid", "attitude", "total"]:
        m = re.search(rf"{prefix}_average_ns (\d+)", combined)
        if m:
            result[f"{prefix}_avg_us"] = round(int(m.group(1)) / 1000, 1)
        m = re.search(rf"{prefix}_95%_ns (\d+)", combined)
        if m:
            result[f"{prefix}_p95_us"] = round(int(m.group(1)) / 1000, 1)

    # Stars identified per match, from the star-ID comparator (images with any identified stars)
    stars_per_match = [s["correct"] + s["incorrect"] for s in starid_by_image.values()
                       if s["correct"] + s["incorrect"] > 0]
    result["avg_stars_per_match"] = (round(sum(stars_per_match) / len(stars_per_match), 1)
                                     if stars_per_match else 0)

    return result


def merge_results(results_list):
    """Merge results from parallel batches."""
    merged = {
        "num_images": sum(r["num_images"] for r in results_list),
        "images_correct": sum(r["images_correct"] for r in results_list),
        "images_incorrect": sum(r["images_incorrect"] for r in results_list),
        "images_too_few_centroids": sum(r["images_too_few_centroids"] for r in results_list),
    }

    all_centroids_avg = [r["avg_centroids"] * r["num_images"] for r in results_list]
    merged["avg_centroids"] = round(sum(all_centroids_avg) / merged["num_images"], 1)
    merged["min_centroids"] = min(r["min_centroids"] for r in results_list)

    all_errors = []
    all_failures = []
    all_per_image = []
    offset = 0
    for r in results_list:
        for p in r.get("per_image", []):
            p["image"] = p["batch_image"] + offset
        offset += r["num_images"]
        all_errors.extend(r.get("attitude_errors_deg", []))
        all_failures.extend(r.get("failures", []))
        all_per_image.extend(r.get("per_image", []))
    merged["attitude_errors_deg"] = all_errors
    merged["failures"] = all_failures
    merged["per_image"] = all_per_image

    if all_errors:
        sorted_errs = sorted(all_errors)
        merged["attitude_error_mean_deg"] = round(sum(all_errors) / len(all_errors), 6)
        merged["attitude_error_median_deg"] = round(sorted_errs[len(sorted_errs) // 2], 6)
        merged["attitude_error_max_deg"] = round(max(all_errors), 6)
        merged["attitude_error_p95_deg"] = round(sorted_errs[int(0.95 * len(sorted_errs))], 6)

    # Average timing across batches (weighted)
    for prefix in ["centroiding", "starid", "attitude", "total"]:
        key_avg = f"{prefix}_avg_us"
        key_p95 = f"{prefix}_p95_us"
        vals_avg = [r[key_avg] for r in results_list if key_avg in r]
        vals_p95 = [r[key_p95] for r in results_list if key_p95 in r]
        if vals_avg:
            merged[key_avg] = round(sum(vals_avg) / len(vals_avg), 1)
        if vals_p95:
            merged[key_p95] = round(max(vals_p95), 1)

    total_matches = sum(
        (r["images_correct"] + r["images_incorrect"]) * r.get("avg_stars_per_match", 0)
        for r in results_list
    )
    total_match_images = sum(r["images_correct"] + r["images_incorrect"] for r in results_list)
    merged["avg_stars_per_match"] = round(total_matches / total_match_images, 1) if total_match_images else 0

    return merged


def format_report(label, config, db_stats, bench, total_wall_time):
    n = bench["num_images"]
    correct = bench["images_correct"]
    incorrect = bench["images_incorrect"]
    too_few = bench["images_too_few_centroids"]
    identifiable = n - too_few
    no_id = identifiable - correct - incorrect

    lines = []
    lines.append(f"{'='*60}")
    lines.append(f"  TETRA BENCHMARK: {label}")
    lines.append(f"  Total wall time: {total_wall_time:.1f}s")
    lines.append(f"  Config: fov={config['fov']}, max-angle={config['max_angle']}, "
                 f"n={n}, parallel={config.get('parallel', 1)}")
    lines.append(f"{'='*60}")

    lines.append(f"\n--- Memory ---")
    lines.append(f"  Database size:          {db_stats['db_size_mb']} MB")
    lines.append(f"  Pattern count:          {db_stats['num_patterns']:,}")

    lines.append(f"\n--- Latency ---")
    lines.append(f"  DB generation:          {db_stats['db_generation_time_s']:.3f} s")
    lines.append(f"  Query wall time:        {bench.get('query_wall_time_s', 0):.3f} s ({n} images)")
    lines.append(f"")
    lines.append(f"  Per-image breakdown (avg / p95):")
    for prefix, label_str in [("centroiding", "Centroiding"), ("starid", "Star-ID"),
                               ("attitude", "Attitude est."), ("total", "Total pipeline")]:
        avg = bench.get(f"{prefix}_avg_us")
        p95 = bench.get(f"{prefix}_p95_us")
        if avg is not None:
            unit = "us" if avg < 10000 else "ms"
            div = 1 if avg < 10000 else 1000
            avg_s, p95_s = avg / div, (p95 or 0) / div
            lines.append(f"    {label_str+':':<22}{avg_s:>8.1f} {unit}  / {p95_s:>8.1f} {unit}")

    lines.append(f"\n--- Star-ID ({n} images, {too_few} had <4 centroids) ---")
    lines.append(f"  Cannot identify:        {too_few}/{n} (insufficient centroids)")
    if identifiable:
        lines.append(f"  Correct:                {correct}/{identifiable} identifiable ({correct/identifiable*100:.1f}%)")
    lines.append(f"  Incorrect IDs:          {incorrect}")
    lines.append(f"  Failed to match:        {no_id}")
    lines.append(f"  Avg centroids/image:    {bench['avg_centroids']}")
    lines.append(f"  Avg stars/match:        {bench['avg_stars_per_match']}")

    errors = bench.get("attitude_errors_deg", [])
    if errors:
        lines.append(f"\n--- Attitude Error ({len(errors)} images with attitude) ---")
        lines.append(f"  Median:                 {bench['attitude_error_median_deg']:.6f} deg")
        lines.append(f"  Mean:                   {bench['attitude_error_mean_deg']:.6f} deg")
        lines.append(f"  P95:                    {bench['attitude_error_p95_deg']:.6f} deg")
        lines.append(f"  Max:                    {bench['attitude_error_max_deg']:.6f} deg")

    failures = bench.get("failures", [])
    if failures:
        lines.append(f"\n--- Failures ({len(failures)} images) ---")
        for f in sorted(failures, key=lambda x: -(x["error_deg"] or 0))[:25]:
            e = f["expected"]
            sid = f.get("starid", {})
            sid_str = f" ids={sid['correct']}ok/{sid['incorrect']}wrong" if sid.get("incorrect") else ""
            if f["error_deg"] is not None:
                a = f["actual"]
                lines.append(f"  image {f['image']:>4d}: {f['error_deg']:>9.4f} deg  "
                             f"expected=({e['ra']:.1f}, {e['de']:.1f}, {e['roll']:.1f})  "
                             f"actual=({a['ra']:.1f}, {a['de']:.1f}, {a['roll']:.1f})  "
                             f"[{f['reason']}]{sid_str}")
            else:
                lines.append(f"  image {f['image']:>4d}:            --  "
                             f"expected=({e['ra']:.1f}, {e['de']:.1f}, {e['roll']:.1f})  "
                             f"[{f['reason']}]{sid_str}")

    lines.append(f"\n{'='*60}")
    return "\n".join(lines)


def run_batch_worker(args):
    """Worker for parallel execution."""
    db_path, fov, batch_size, seed, algo, extra_flags = args
    output = run_batch(db_path, fov, batch_size, seed=seed, star_id_algo=algo, extra_flags=extra_flags)
    return parse_results(output, batch_size, seed)


def main():
    parser = argparse.ArgumentParser(description="Tetra benchmark")
    parser.add_argument("--num-images", type=int, default=100)
    parser.add_argument("--fov", type=float, default=20)
    parser.add_argument("--max-angle", type=float, default=None,
                        help="Tetra database max angle. Defaults to fov + 5")
    parser.add_argument("--algo", type=str, default="tetra",
                        choices=["tetra", "py", "both"],
                        help="Star-id algorithm: tetra, py (pyramid), or both")
    parser.add_argument("--label", type=str, default="baseline")
    parser.add_argument("--parallel", type=int, default=4,
                        help="Number of parallel LOST processes")
    parser.add_argument("--tetra-pairing", type=str, default=None,
                        choices=["triad", "centroid-radius"],
                        help="Tetra star pairing mode")
    parser.add_argument("--db-flags", type=str, default="",
                        help="Extra flags for Tetra database generation, e.g. '--min-mag 5'")
    parser.add_argument("--output", type=str, default=None,
                        help="Save results as JSON to this path")
    args = parser.parse_args()

    if args.max_angle is None:
        args.max_angle = args.fov + 5

    wall_t0 = time.monotonic()
    db_path = DB_PATH_TEMPLATE.format(label=args.label)

    config = {
        "num_images": args.num_images,
        "fov": args.fov,
        "max_angle": args.max_angle,
        "parallel": args.parallel,
    }

    algos = ["tetra", "py"] if args.algo == "both" else [args.algo]
    all_reports = []

    extra_flags = ""
    if args.tetra_pairing:
        extra_flags += f" --tetra-pairing {args.tetra_pairing}"
        config["tetra_pairing"] = args.tetra_pairing
    if args.db_flags:
        config["db_flags"] = args.db_flags

    for algo in algos:
        algo_label = f"{args.label}-{algo}" if len(algos) > 1 else args.label
        algo_name = "Tetra" if algo == "tetra" else "Pyramid"
        algo_db_path = DB_PATH_TEMPLATE.format(label=algo_label)

        # Step 1: Generate database
        print(f"\n[{algo_name}] Generating database...")
        db_stats = generate_database(algo, args.max_angle, algo_db_path, args.db_flags)
        print(f"  Done: {db_stats['db_size_mb']} MB, {db_stats['db_generation_time_s']:.1f}s"
              + (f", {db_stats['num_patterns']:,} patterns" if db_stats['num_patterns'] else ""))

        # Step 2: Run images (possibly in parallel batches)
        num_workers = args.parallel
        batch_size = args.num_images // num_workers
        remainder = args.num_images % num_workers

        batches = []
        for i in range(num_workers):
            n = batch_size + (1 if i < remainder else 0)
            batches.append((algo_db_path, args.fov, n, 42 + i, algo, extra_flags))

        print(f"[{algo_name}] Running {args.num_images} images ({num_workers} batch(es))...")
        t0 = time.monotonic()

        if num_workers == 1:
            output = run_batch(algo_db_path, args.fov, args.num_images, seed=42, star_id_algo=algo, extra_flags=extra_flags)
            bench = parse_results(output, args.num_images, 42)
        else:
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures = [executor.submit(run_batch_worker, b) for b in batches]
                results = [f.result() for f in futures]
            bench = merge_results(results)

        query_time = time.monotonic() - t0
        bench["query_wall_time_s"] = round(query_time, 3)

        print(f"  Done: {bench['images_correct']}/{args.num_images} correct, "
              f"{query_time:.1f}s wall time")

        algo_wall = time.monotonic() - wall_t0

        # Print and save report
        report = format_report(f"{args.label} ({algo_name})", config, db_stats, bench, algo_wall)
        print(f"\n{report}\n")
        all_reports.append(report)

        # Save full results per algo
        output_path = os.path.join(BENCH_DIR, f"benchmark-{algo_label}.json")
        full_results = {
            "label": algo_label,
            "algo": algo,
            "config": config,
            "database": db_stats,
            "benchmark": {k: v for k, v in bench.items() if k != "attitude_errors_deg"},
            "per_image": bench.get("per_image", []),
        }
        with open(output_path, "w") as f:
            json.dump(full_results, f, indent=2)
        print(f"Full results saved to {output_path}")

        # Genuine failures only (wrong IDs/attitude or failed to match), with exact repro commands.
        # Image batch_image of a run with seed S is regenerated by --generate (batch_image+1)
        # --generate-seed S; it is the last image of that run.
        genuine = []
        for p in bench.get("per_image", []):
            if p["reason"] in ("correct", "too_few_centroids"):
                continue
            genuine.append({**p, "algo": algo, "repro": pipeline_cmd(
                algo_db_path, args.fov, p["batch_image"] + 1, p["seed"], algo, extra_flags)})
        genuine_path = os.path.join(BENCH_DIR, f"failures-{algo_label}.json")
        with open(genuine_path, "w") as f:
            json.dump(genuine, f, indent=2)
        print(f"{len(genuine)} genuine failures saved to {genuine_path}")

        # Accumulate new failures into shared failure-cases.json
        failure_path = os.path.join(BENCH_DIR, "failure-cases.json")
        existing_failures = []
        if os.path.exists(failure_path):
            with open(failure_path) as f:
                existing_failures = json.load(f)

        new_failures = []
        for p in bench.get("per_image", []):
            if p["reason"] == "correct":
                continue
            entry = {**p, "algo": algo, "source": args.label}
            if extra_flags:
                entry["flags"] = extra_flags.strip()
            new_failures.append(entry)

        if new_failures:
            # Deduplicate by (algo, expected attitude, reason)
            existing_keys = set()
            for ef in existing_failures:
                e = ef.get("expected", {})
                key = (ef.get("algo"), round(e.get("ra",0),1), round(e.get("de",0),1),
                       round(e.get("roll",0),1), ef.get("reason"))
                existing_keys.add(key)

            added = 0
            for nf in new_failures:
                e = nf.get("expected", {})
                key = (nf.get("algo"), round(e.get("ra",0),1), round(e.get("de",0),1),
                       round(e.get("roll",0),1), nf.get("reason"))
                if key not in existing_keys:
                    existing_failures.append(nf)
                    existing_keys.add(key)
                    added += 1

            with open(failure_path, "w") as f:
                json.dump(existing_failures, f, indent=2)
            if added:
                print(f"Added {added} new failure cases to {failure_path} (total: {len(existing_failures)})")

    # Save combined report
    os.makedirs(BENCH_DIR, exist_ok=True)
    report_path = os.path.join(BENCH_DIR, f"benchmark-{args.label}.txt")
    with open(report_path, "w") as f:
        f.write("\n\n".join(all_reports) + "\n")
    print(f"\nReport saved to {report_path}")


if __name__ == "__main__":
    main()
