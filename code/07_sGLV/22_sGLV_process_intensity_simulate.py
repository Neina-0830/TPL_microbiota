#!/usr/bin/env python3
"""
Generate revised top-10 process-intensity sGLV simulations.

Each TPL formulation uses its own ten lowest-loss calibrated parameter sets.
For each retained set and condition, one focal process is varied over
z = 0, 0.25, 0.50, 0.75, 1 while all other components remain calibrated.

Examples
--------
python code/07_sGLV/22_sGLV_process_intensity_simulate.py \
  --tpl tTTPL --group all --condition both

python code/07_sGLV/22_sGLV_process_intensity_simulate.py \
  --tpl pTTPL --group Bacteria --process immigration noise
"""
from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

from sGLV_core import (
    CONDITIONS, INTENSITIES, PROCESSES, SEEDS, load_condition_inputs,
    load_group_config, process_filename, process_intensity_parameters,
    project_root, work_root, simulate_sglv, validate_group_names, write_simulation
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tpl", choices=["tTTPL", "pTTPL"], required=True)
    p.add_argument("--group", nargs="+", default=["all"])
    p.add_argument("--condition", choices=["warming", "unwarming", "both"], default="both")
    p.add_argument("--process", nargs="+", default=["all"])
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--steps-per-snapshot", type=int, default=13)
    p.add_argument("--n-snapshots", type=int, default=12)
    return p.parse_args()


def top10_file(work: Path, tpl: str, group: str) -> Path:
    return work / "Calibration" / tpl / f"{group}_top10_parameters.csv"


def main():
    args = parse_args()
    root = project_root()
    work = work_root(root)
    configs = load_group_config()
    groups = validate_group_names(args.group, configs)
    conditions = CONDITIONS if args.condition == "both" else (args.condition,)

    if args.process == ["all"]:
        processes = PROCESSES
    else:
        bad = [p for p in args.process if p not in PROCESSES]
        if bad:
            raise ValueError(f"Unknown process(es): {bad}")
        processes = tuple(args.process)

    for group in groups:
        cfg = configs[group]
        pfile = top10_file(work, args.tpl, group)
        if not pfile.exists():
            raise FileNotFoundError(
                f"Missing top-10 calibration file: {pfile}\n"
                "Run 21_sGLV_select_top10.R first."
            )
        top = pd.read_csv(pfile)
        need = {"Rank", "Immigration", "Interaction", "Growth", "Noise"}
        miss = need - set(top.columns)
        if miss:
            raise ValueError(f"{pfile}: missing columns {sorted(miss)}")
        if len(top) != 10:
            raise ValueError(f"{pfile}: expected 10 retained parameter sets, found {len(top)}.")

        for condition in conditions:
            current_a = None
            inp = None
            for _, row in top.sort_values("Rank").iterrows():
                rank = int(row["Rank"])
                d = float(row["Immigration"])
                a = float(row["Interaction"])
                g = float(row["Growth"])
                noise = float(row["Noise"])
                if current_a is None or a != current_a:
                    inp = load_condition_inputs(root, cfg, condition, a)
                    current_a = a

                for process in processes:
                    outdir = (
                        work / "get_timeseries" / group / args.tpl /
                        "process_intensity_top10" / process / condition
                    )
                    for z in INTENSITIES:
                        for seed in SEEDS:
                            path = outdir / process_filename(
                                rank, d, a, g, seed, noise, process, z
                            )
                            if path.exists() and not args.overwrite:
                                continue
                            if args.dry_run:
                                print(path)
                                continue
                            params = process_intensity_parameters(
                                inputs=inp, cfg=cfg, condition=condition,
                                immigration=d, interaction=a, growth=g, noise=noise,
                                process=process, intensity=z,
                            )
                            out = simulate_sglv(
                                params=params,
                                taxa=inp["taxa"],
                                seed=seed,
                                clip_max=cfg["resample_depth"],
                                steps_per_snapshot=args.steps_per_snapshot,
                                n_snapshots=args.n_snapshots,
                            )
                            write_simulation(out, path, overwrite=args.overwrite)
            print(f"Finished {args.tpl} | {group} | {condition}")


if __name__ == "__main__":
    main()
