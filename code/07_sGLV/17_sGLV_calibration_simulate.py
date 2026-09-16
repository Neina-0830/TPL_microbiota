#!/usr/bin/env python3
"""
Generate the full calibration-grid sGLV simulations for all five soil-biotic groups.

Examples
--------
# List the workload without running:
python code/07_sGLV/17_sGLV_calibration_simulate.py --group all --condition both --dry-run

# Run Archaea warming:
python code/07_sGLV/17_sGLV_calibration_simulate.py --group Archaea --condition warming

# HPC chunking (0-based slice of the parameter/seed grid):
python code/07_sGLV/17_sGLV_calibration_simulate.py \
  --group Bacteria --condition warming --start-index 0 --stop-index 500

All external parameter inputs are read from ./data/sGLV by default.
"""
from __future__ import annotations
import argparse
from pathlib import Path

from sGLV_core import (
    CONDITIONS, calibration_filename, calibration_grid, calibrated_parameters,
    load_condition_inputs, load_group_config, project_root, work_root, simulate_sglv,
    validate_group_names, write_simulation
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--group", nargs="+", default=["all"])
    p.add_argument("--condition", choices=["warming", "unwarming", "both"], default="both")
    p.add_argument("--start-index", type=int, default=0,
                   help="0-based inclusive grid index for HPC chunking")
    p.add_argument("--stop-index", type=int, default=None,
                   help="0-based exclusive grid index for HPC chunking")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--steps-per-snapshot", type=int, default=13)
    p.add_argument("--n-snapshots", type=int, default=12)
    return p.parse_args()


def main():
    args = parse_args()
    root = project_root()
    out_root = work_root(root)
    configs = load_group_config()
    groups = validate_group_names(args.group, configs)
    conditions = CONDITIONS if args.condition == "both" else (args.condition,)

    for group in groups:
        cfg = configs[group]
        grid = calibration_grid(cfg)
        start = max(0, args.start_index)
        stop = len(grid) if args.stop_index is None else min(len(grid), args.stop_index)
        chunk = grid[start:stop]

        for condition in conditions:
            outdir = out_root / "get_timeseries" / group / "raw" / condition
            print(f"[{group} | {condition}] total grid={len(grid)}, this run={len(chunk)} "
                  f"(indices {start}:{stop})")
            if args.dry_run:
                continue

            # Hold at most one potentially very large interaction matrix in memory.
            current_a = None
            inp = None
            for k, (d, a, g, seed, noise) in enumerate(chunk, start=start):
                if current_a is None or a != current_a:
                    inp = load_condition_inputs(root, cfg, condition, a)
                    current_a = a
                params = calibrated_parameters(inp, cfg, condition, d, a, g, noise)
                out = simulate_sglv(
                    params=params,
                    taxa=inp["taxa"],
                    seed=seed,
                    clip_max=cfg["resample_depth"],
                    steps_per_snapshot=args.steps_per_snapshot,
                    n_snapshots=args.n_snapshots,
                )
                path = outdir / calibration_filename(d, a, g, seed, noise)
                wrote = write_simulation(out, path, overwrite=args.overwrite)
                if (k - start + 1) % 50 == 0 or k == stop - 1:
                    print(f"  {k+1}/{stop}: {'wrote' if wrote else 'exists'} {path.name}")


if __name__ == "__main__":
    main()
