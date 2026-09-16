#!/usr/bin/env python3
"""Preflight check for external sGLV parameter inputs under ./data/sGLV.

This script does not create or infer model parameterizations.  It verifies that
all condition-specific empirical vectors and interaction matrices needed by the
calibration and process-intensity simulations are present before a large HPC run.

Run from the repository root:
    python code/07_sGLV/16_check_sGLV_inputs.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

from sGLV_core import (
    CONDITIONS,
    data_root,
    initial_abundance_file,
    immigration_file,
    growth_file,
    interaction_candidates,
    load_group_config,
    project_root,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--group", nargs="+", default=["all"],
                   help="Group(s) to check; default: all")
    p.add_argument("--quiet", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = project_root()
    cfgs = load_group_config()
    groups = list(cfgs) if args.group == ["all"] else args.group
    unknown = [g for g in groups if g not in cfgs]
    if unknown:
        raise ValueError(f"Unknown group(s): {unknown}")

    droot = data_root(root)
    missing: list[tuple[str, str]] = []
    found = 0

    for group in groups:
        cfg = cfgs[group]
        for condition in CONDITIONS:
            basics = [
                ("initial abundance", initial_abundance_file(root, cfg, condition)),
                ("relative immigration", immigration_file(root, cfg, condition)),
                ("relative growth", growth_file(root, cfg, condition)),
            ]
            for label, path in basics:
                if path.exists():
                    found += 1
                else:
                    missing.append((f"{group} | {condition} | {label}", str(path)))

            for strength in cfg["interaction_grid"]:
                candidates = interaction_candidates(root, cfg, condition, strength)
                if any(p.exists() for p in candidates):
                    found += 1
                else:
                    missing.append((
                        f"{group} | {condition} | interaction={strength}",
                        " OR ".join(str(p) for p in candidates),
                    ))

    if not args.quiet:
        print(f"Repository root: {root}")
        print(f"sGLV input root: {droot}")
        print(f"Found required input items: {found}")
        print(f"Missing required input items: {len(missing)}")

    if missing:
        print("\nMissing sGLV model-derived inputs:")
        for label, path in missing:
            print(f"- {label}\n    {path}")
        print(
            "\nThese files encode the empirical immigration, growth and interaction "
            "structures used to parameterize the sGLV model. They are not reconstructed "
            "automatically from the taxon tables because doing so would change the "
            "scientific parameterization. Place the original files under ./data/sGLV "
            "using the layout documented in data/sGLV/README.md."
        )
        return 1

    print("\nAll requested sGLV external parameter inputs are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
