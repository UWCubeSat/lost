#!/usr/bin/env python3
"""
Noise sweep benchmark for LOST star tracker.

Runs Tetra and Pyramid across independent sweeps of read noise, false stars,
and centroid perturbation. Produces per-config JSON results and a summary
comparison table.

Usage:
    python3 noise_benchmark.py [--num-images 1000] [--fov 20] [--label noise-sweep]
                               [--max-concurrent 4] [--inner-parallel 2]
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

from benchmark import (
    BENCH_DIR,
    DB_PATH_TEMPLATE,
    generate_database,
    merge_results,
    parse_results,
    run_batch,
    run_batch_worker,
)

ALGOS = ["tetra", "py"]
ALGO_NAMES = {"tetra": "Tetra", "py": "Pyramid"}

# Sweep definitions: each is (parameter_name, cli_flag, levels)
# The first level in each sweep is the baseline value.
SWEEPS = {
    "read_noise": {
        "flag": "--generate-read-noise-stddev",
        "levels": [0, 0.01, 0.05, 0.10, 0.15],
        "baseline": 0.05,
        "unit": "",
    },
    "false_stars": {
        "flag": "--generate-false-stars",
        "levels": [0, 100, 500],
        "baseline": 0,
        "unit": "",
    },
    "perturb_centroids": {
        "flag": "--generate-perturb-centroids",
        "levels": [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
        "baseline": 0,
        "unit": "px",
    },
}

# Baseline generate flags (shared across all runs unless overridden by sweep)
BASELINE_FLAGS = {
    "--generate-read-noise-stddev": "0.05",
    "--generate-false-stars": "0",
    "--generate-perturb-centroids": "0",
}


def build_generate_flags(sweep_name, level):
    """Build the extra generate flags for a sweep config."""
    flags = dict(BASELINE_FLAGS)
    flags[SWEEPS[sweep_name]["flag"]] = str(level)
    # For the zero-noise case, also disable shot noise and dark current
    if sweep_name == "read_noise" and level == 0:
        flags["--generate-shot-noise"] = "false"
        flags["--generate-dark-current"] = "0"
    return " ".join(f"{k} {v}" for k, v in flags.items())


def config_label(base_label, sweep_name, level):
    return f"{base_label}-{sweep_name}-{level}"


def run_single_config(args_tuple):
    """Run a single (algo, sweep, level) config. Intended for process pool."""
    (algo, sweep_name, level, num_images, fov, max_angle, inner_parallel,
     base_label, db_flags) = args_tuple

    label = config_label(base_label, sweep_name, level)
    algo_label = f"{label}-{algo}"
    db_path = DB_PATH_TEMPLATE.format(label=algo_label)
    extra_flags = build_generate_flags(sweep_name, level)

    db_stats = generate_database(algo, max_angle, db_path, db_flags)

    num_workers = inner_parallel
    batch_size = num_images // num_workers
    remainder = num_images % num_workers

    if num_workers == 1:
        output = run_batch(db_path, fov, num_images, seed=42,
                           star_id_algo=algo, extra_flags=extra_flags)
        bench = parse_results(output, num_images, 42)
    else:
        batches = []
        for i in range(num_workers):
            n = batch_size + (1 if i < remainder else 0)
            batches.append((db_path, fov, n, 42 + i, algo, extra_flags))

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(run_batch_worker, b) for b in batches]
            results = [f.result() for f in futures]
        bench = merge_results(results)

    return {
        "algo": algo,
        "algo_name": ALGO_NAMES[algo],
        "sweep": sweep_name,
        "level": level,
        "label": algo_label,
        "db_stats": db_stats,
        "bench": bench,
    }


def format_comparison_table(results_by_sweep):
    """Format a summary comparison table across all sweeps."""
    lines = []

    for sweep_name, sweep_results in results_by_sweep.items():
        sweep_def = SWEEPS[sweep_name]
        unit = f" {sweep_def['unit']}" if sweep_def["unit"] else ""
        lines.append(f"\n{'='*100}")
        lines.append(f"  SWEEP: {sweep_name} ({sweep_def['flag']})")
        lines.append(f"{'='*100}")

        header = (f"{'Level':>10} | {'Algo':<8} | {'Correct':>9} | {'Incorr':>6} | "
                  f"{'Centroid':>10} | {'StarID':>10} | {'Attitude':>10} | {'Total':>10} | "
                  f"{'Att err med':>12} | {'Att err p95':>12}")
        lines.append(header)
        lines.append("-" * len(header))

        for r in sorted(sweep_results, key=lambda x: (x["level"], x["algo"])):
            b = r["bench"]
            n = b["num_images"]
            too_few = b["images_too_few_centroids"]
            identifiable = n - too_few
            correct = b["images_correct"]
            incorrect = b["images_incorrect"]

            def fmt_us(key):
                v = b.get(key)
                if v is None:
                    return "--"
                if v < 10000:
                    return f"{v:.0f} us"
                return f"{v/1000:.1f} ms"

            att_med = b.get("attitude_error_median_deg")
            att_p95 = b.get("attitude_error_p95_deg")

            lines.append(
                f"{str(r['level'])+unit:>10} | {r['algo_name']:<8} | "
                f"{correct:>4}/{identifiable:<4} | {incorrect:>6} | "
                f"{fmt_us('centroiding_avg_us'):>10} | {fmt_us('starid_avg_us'):>10} | "
                f"{fmt_us('attitude_avg_us'):>10} | {fmt_us('total_avg_us'):>10} | "
                f"{att_med:>11.4f}°" if att_med is not None else f"{'--':>12}" + " | "
                f"{att_p95:>11.4f}°" if att_p95 is not None else f"{'--':>12}"
            )

        lines.append("")

    return "\n".join(lines)


def format_comparison_table(results_by_sweep):
    """Format a summary comparison table across all sweeps."""
    lines = []

    for sweep_name, sweep_results in results_by_sweep.items():
        sweep_def = SWEEPS[sweep_name]
        unit = f" {sweep_def['unit']}" if sweep_def["unit"] else ""
        lines.append(f"\n{'='*110}")
        lines.append(f"  SWEEP: {sweep_name} ({sweep_def['flag']})")
        lines.append(f"{'='*110}")

        header = (f"  {'Level':>10} | {'Algo':<8} | {'Correct':>9} | {'Incorr':>6} | "
                  f"{'Centroid':>10} | {'StarID':>10} | {'Attitude':>10} | {'Total':>10} | "
                  f"{'Att med':>10} | {'Att p95':>10}")
        lines.append(header)
        lines.append("  " + "-" * 106)

        for r in sorted(sweep_results, key=lambda x: (x["level"], x["algo"])):
            b = r["bench"]
            n = b["num_images"]
            too_few = b["images_too_few_centroids"]
            identifiable = n - too_few
            correct = b["images_correct"]
            incorrect = b["images_incorrect"]

            def fmt_us(key):
                v = b.get(key)
                if v is None:
                    return "--"
                if v < 10000:
                    return f"{v:.0f} us"
                return f"{v/1000:.1f} ms"

            att_med = b.get("attitude_error_median_deg")
            att_p95 = b.get("attitude_error_p95_deg")
            att_med_s = f"{att_med:.4f} deg" if att_med is not None else "--"
            att_p95_s = f"{att_p95:.4f} deg" if att_p95 is not None else "--"

            level_s = f"{r['level']}{unit}"
            lines.append(
                f"  {level_s:>10} | {r['algo_name']:<8} | "
                f"{correct:>4}/{identifiable:<4} | {incorrect:>6} | "
                f"{fmt_us('centroiding_avg_us'):>10} | {fmt_us('starid_avg_us'):>10} | "
                f"{fmt_us('attitude_avg_us'):>10} | {fmt_us('total_avg_us'):>10} | "
                f"{att_med_s:>10} | {att_p95_s:>10}"
            )

        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="LOST noise sweep benchmark")
    parser.add_argument("--num-images", type=int, default=1000)
    parser.add_argument("--fov", type=float, default=20)
    parser.add_argument("--max-angle", type=float, default=None,
                        help="Tetra database max angle. Defaults to fov + 5")
    parser.add_argument("--label", type=str, default="noise-sweep")
    parser.add_argument("--max-concurrent", type=int, default=4,
                        help="Max concurrent benchmark configs to run")
    parser.add_argument("--inner-parallel", type=int, default=2,
                        help="Parallel batches within each config")
    parser.add_argument("--db-flags", type=str, default="",
                        help="Extra flags for database generation")
    parser.add_argument("--sweeps", type=str, default="all",
                        help="Comma-separated sweep names to run, or 'all'")
    args = parser.parse_args()

    if args.max_angle is None:
        args.max_angle = args.fov + 5

    sweep_names = (list(SWEEPS.keys()) if args.sweeps == "all"
                   else args.sweeps.split(","))

    # Build the list of all (algo, sweep, level) configs, deduplicating baseline
    configs = []
    seen = set()
    for sweep_name in sweep_names:
        for level in SWEEPS[sweep_name]["levels"]:
            for algo in ALGOS:
                config_key = (algo, sweep_name, level)
                if config_key not in seen:
                    seen.add(config_key)
                    configs.append((
                        algo, sweep_name, level, args.num_images, args.fov,
                        args.max_angle, args.inner_parallel, args.label,
                        args.db_flags,
                    ))

    total = len(configs)
    print(f"Noise sweep: {total} configs ({len(sweep_names)} sweeps x 2 algos)")
    print(f"  {args.num_images} images/config, fov={args.fov}, "
          f"max-concurrent={args.max_concurrent}, inner-parallel={args.inner_parallel}")
    print(f"  Sweeps: {', '.join(sweep_names)}")
    print()

    wall_t0 = time.monotonic()
    all_results = []
    completed = 0

    with ProcessPoolExecutor(max_workers=args.max_concurrent) as executor:
        future_to_config = {}
        for cfg in configs:
            algo, sweep_name, level = cfg[0], cfg[1], cfg[2]
            f = executor.submit(run_single_config, cfg)
            future_to_config[f] = (algo, sweep_name, level)

        for future in as_completed(future_to_config):
            algo, sweep_name, level = future_to_config[future]
            completed += 1
            try:
                result = future.result()
                all_results.append(result)
                b = result["bench"]
                print(f"  [{completed}/{total}] {ALGO_NAMES[algo]:>7} | "
                      f"{sweep_name}={level}: "
                      f"{b['images_correct']}/{b['num_images']} correct, "
                      f"{b['images_incorrect']} incorrect")
            except Exception as e:
                print(f"  [{completed}/{total}] {ALGO_NAMES[algo]:>7} | "
                      f"{sweep_name}={level}: FAILED: {e}", file=sys.stderr)

    wall_time = time.monotonic() - wall_t0
    print(f"\nAll configs done in {wall_time:.1f}s")

    # Group results by sweep
    results_by_sweep = {}
    for r in all_results:
        results_by_sweep.setdefault(r["sweep"], []).append(r)

    # Print comparison table
    table = format_comparison_table(results_by_sweep)
    print(table)

    # Save per-config JSONs
    os.makedirs(BENCH_DIR, exist_ok=True)
    for r in all_results:
        b = r["bench"]
        output = {
            "label": r["label"],
            "algo": r["algo"],
            "sweep": r["sweep"],
            "level": r["level"],
            "config": {
                "num_images": b["num_images"],
                "fov": args.fov,
                "max_angle": args.max_angle,
            },
            "database": r["db_stats"],
            "benchmark": {k: v for k, v in b.items()
                          if k not in ("attitude_errors_deg", "per_image", "failures")},
            "attitude_errors_deg": b.get("attitude_errors_deg", []),
        }
        path = os.path.join(BENCH_DIR, f"benchmark-{r['label']}.json")
        with open(path, "w") as f:
            json.dump(output, f, indent=2)

    # Save combined summary
    summary = {
        "label": args.label,
        "wall_time_s": round(wall_time, 1),
        "config": {
            "num_images": args.num_images,
            "fov": args.fov,
            "max_angle": args.max_angle,
            "max_concurrent": args.max_concurrent,
            "inner_parallel": args.inner_parallel,
        },
        "sweeps": {},
    }
    for sweep_name, sweep_results in results_by_sweep.items():
        summary["sweeps"][sweep_name] = []
        for r in sorted(sweep_results, key=lambda x: (x["level"], x["algo"])):
            b = r["bench"]
            n = b["num_images"]
            too_few = b["images_too_few_centroids"]
            identifiable = n - too_few
            summary["sweeps"][sweep_name].append({
                "level": r["level"],
                "algo": r["algo"],
                "correct": b["images_correct"],
                "identifiable": identifiable,
                "incorrect": b["images_incorrect"],
                "centroiding_avg_us": b.get("centroiding_avg_us"),
                "starid_avg_us": b.get("starid_avg_us"),
                "attitude_avg_us": b.get("attitude_avg_us"),
                "total_avg_us": b.get("total_avg_us"),
                "centroiding_p95_us": b.get("centroiding_p95_us"),
                "starid_p95_us": b.get("starid_p95_us"),
                "attitude_p95_us": b.get("attitude_p95_us"),
                "total_p95_us": b.get("total_p95_us"),
                "attitude_error_median_deg": b.get("attitude_error_median_deg"),
                "attitude_error_mean_deg": b.get("attitude_error_mean_deg"),
                "attitude_error_p95_deg": b.get("attitude_error_p95_deg"),
                "attitude_error_max_deg": b.get("attitude_error_max_deg"),
            })

    summary_path = os.path.join(BENCH_DIR, f"benchmark-{args.label}-summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to {summary_path}")

    # Save the table as text
    table_path = os.path.join(BENCH_DIR, f"benchmark-{args.label}.txt")
    with open(table_path, "w") as f:
        f.write(f"LOST Noise Sweep Benchmark: {args.label}\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Config: {args.num_images} images, fov={args.fov}, "
                f"wall time={wall_time:.1f}s\n")
        f.write(table + "\n")
    print(f"Report saved to {table_path}")


if __name__ == "__main__":
    main()
