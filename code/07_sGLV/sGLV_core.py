#!/usr/bin/env python3
"""
Shared sGLV utilities for the TPL process-sensitivity workflow.

This refactor keeps the ecological parameterization of the supplied all-group
simulation scripts but removes machine-specific paths and the dependency on the
project-local ``noise_parameters.py`` module.  Only the GLV + multiplicative
linear Langevin implementation used by the manuscript workflow is retained.

Repository layout
-----------------
Run from the repository root. Primary sGLV parameter inputs are read from
``./data/sGLV`` (override with ``TPL_SGLV_DATA_DIR``), while generated
simulation/calibration products are written under ``./results/17_sGLV``.

The historical output prefix ``SLM`` is retained only for filename compatibility.
The model referred to in the manuscript is the stochastic generalized
Lotka--Volterra (sGLV) model.
"""
from __future__ import annotations

import csv
import itertools
import os
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CONFIG_FILE = HERE / "sGLV_group_config.csv"

CONDITIONS = ("warming", "unwarming")
PROCESSES = ("immigration", "growth", "interaction", "noise")
SEEDS = (0, 1, 2)
INTENSITIES = (0.0, 0.25, 0.50, 0.75, 1.0)


def project_root() -> Path:
    """Return the repository root."""
    return Path(os.environ.get("TPL_PROJECT_ROOT", os.getcwd())).expanduser().resolve()


def data_root(root: Path | None = None) -> Path:
    """Directory containing external sGLV parameter inputs."""
    root = project_root() if root is None else Path(root)
    return Path(os.environ.get("TPL_SGLV_DATA_DIR", root / "data" / "sGLV")).expanduser().resolve()


def work_root(root: Path | None = None) -> Path:
    """Directory containing generated sGLV simulations and calibration products."""
    root = project_root() if root is None else Path(root)
    return Path(os.environ.get("TPL_SGLV_WORK_DIR", root / "results" / "17_sGLV")).expanduser().resolve()


def _parse_grid(x: str) -> List[float]:
    return [float(v) for v in str(x).split(";") if str(v).strip()]


def load_group_config(config_file: Path = CONFIG_FILE) -> Dict[str, dict]:
    tab = pd.read_csv(config_file)
    ans: Dict[str, dict] = {}
    for _, row in tab.iterrows():
        g = str(row["Group"])
        ans[g] = {
            "group": g,
            "resample_depth": int(row["ResampleDepth"]),
            "immigration_grid": _parse_grid(row["ImmigrationGrid"]),
            "interaction_grid": _parse_grid(row["InteractionGrid"]),
            "growth_grid": _parse_grid(row["GrowthGrid"]),
            "noise_grid": _parse_grid(row["NoiseGrid"]),
            "warmed_drift": float(row["WarmedDrift"]),
            "unwarmed_drift": float(row["UnwarmedDrift"]),
            "warmed_noise_multiplier": float(row["WarmedNoiseMultiplier"]),
        }
    return ans


def numeric_token(x: float) -> str:
    """Stable compact numeric token compatible with the legacy filenames."""
    x = float(x)
    if x == 0:
        return "0"
    # General format reproduces familiar tokens (1e-08, 0.0001, 0.05, 1).
    return f"{x:.12g}"


def intensity_token(z: float) -> str:
    return numeric_token(z).replace(".", "p")


def calibration_filename(d: float, a: float, g: float, seed: int, noise: float) -> str:
    return (
        "Time_series_from_SLM_"
        f"immigration_d{numeric_token(d)}_"
        f"interaction_sd{numeric_token(a)}_"
        f"growth_g{numeric_token(g)}_"
        f"seed{int(seed)}_"
        f"noise{numeric_token(noise)}.csv"
    )


def process_filename(rank: int, d: float, a: float, g: float, seed: int,
                     noise: float, process: str, intensity: float) -> str:
    return (
        f"Time_series_from_SLM_top10rank{int(rank)}_"
        f"immigration_d{numeric_token(d)}_"
        f"interaction_sd{numeric_token(a)}_"
        f"growth_g{numeric_token(g)}_"
        f"seed{int(seed)}_"
        f"noise{numeric_token(noise)}_"
        f"process_{process}_"
        f"intensity{intensity_token(intensity)}.csv"
    )


def initial_abundance_file(root: Path, cfg: Mapping, condition: str) -> Path:
    return data_root(root) / "rawdata" / (
        f"resample{cfg['resample_depth']}_otutab_{cfg['group']}_"
        f"{condition}_2010_mean.csv"
    )


def immigration_file(root: Path, cfg: Mapping, condition: str) -> Path:
    return data_root(root) / "rawdata" / (
        f"otutab_{cfg['group']}_2010to2020_{condition}_relative_immigration.csv"
    )


def growth_file(root: Path, cfg: Mapping, condition: str) -> Path:
    return data_root(root) / "rawdata" / (
        f"otutab_{cfg['group']}_2010to2020_{condition}_copynumber_relative_growth.csv"
    )


def interaction_candidates(root: Path, cfg: Mapping, condition: str, strength: float) -> List[Path]:
    token = numeric_token(strength)
    name = (
        f"otutab_{cfg['group']}_{condition}_mean0_sd{token}_interaction_matrix.csv"
    )
    return [
        data_root(root) / "interaction_matrix" / "data1_data2" / name,
        data_root(root) / "interaction_matrix" / name,
        data_root(root) / "interaction" / name,
    ]


def resolve_existing(candidates: Sequence[Path], label: str) -> Path:
    for p in candidates:
        if p.exists():
            return p
    msg = "\n".join(str(p) for p in candidates)
    raise FileNotFoundError(f"Could not find {label}. Tried:\n{msg}")


def _read_one_numeric_column(path: Path, preferred: str) -> pd.Series:
    dat = pd.read_csv(path, index_col=0)
    if preferred not in dat.columns:
        if dat.shape[1] != 1:
            raise ValueError(
                f"{path}: expected column '{preferred}' or exactly one data column; "
                f"found {list(dat.columns)}"
            )
        col = dat.columns[0]
    else:
        col = preferred
    x = pd.to_numeric(dat[col], errors="raise")
    x.index = x.index.astype(str)
    return x


def load_condition_inputs(root: Path, cfg: Mapping, condition: str,
                          interaction_strength: float) -> dict:
    """
    Load and align the empirical vectors/matrix used by one sGLV condition.

    The taxon order is taken from the initial-abundance file. Immigration,
    growth, and interaction inputs are reordered to the same taxa when labels
    are available.
    """
    if condition not in CONDITIONS:
        raise ValueError(f"Unknown condition: {condition}")

    init_path = initial_abundance_file(root, cfg, condition)
    imm_path = immigration_file(root, cfg, condition)
    gro_path = growth_file(root, cfg, condition)
    int_path = resolve_existing(
        interaction_candidates(root, cfg, condition, interaction_strength),
        f"{cfg['group']} {condition} interaction matrix at {interaction_strength}",
    )

    initial = _read_one_numeric_column(init_path, "MeanAbun")
    immigration = _read_one_numeric_column(imm_path, "Relative_immigration")
    growth_col = f"{condition}_relative_growth"
    growth = _read_one_numeric_column(gro_path, growth_col)

    taxa = initial.index.astype(str)
    if not set(taxa).issubset(set(immigration.index)):
        missing = list(set(taxa) - set(immigration.index))[:10]
        raise ValueError(f"Immigration file is missing taxa, e.g. {missing}")
    if not set(taxa).issubset(set(growth.index)):
        missing = list(set(taxa) - set(growth.index))[:10]
        raise ValueError(f"Growth file is missing taxa, e.g. {missing}")

    immigration = immigration.loc[taxa]
    growth = growth.loc[taxa]

    mat_df = pd.read_csv(int_path, index_col=0)
    # If the matrix has explicit matching row/column taxon labels, align them.
    row_labels = mat_df.index.astype(str)
    col_labels = pd.Index([str(x) for x in mat_df.columns])
    if set(taxa).issubset(set(row_labels)) and set(taxa).issubset(set(col_labels)):
        mat_df.index = row_labels
        mat_df.columns = col_labels
        mat = mat_df.loc[taxa, taxa].to_numpy(dtype=float)
    else:
        mat = mat_df.to_numpy(dtype=float)
        n = len(taxa)
        if mat.size != n * n:
            raise ValueError(
                f"{int_path}: matrix has {mat.shape}, cannot reshape to {n} x {n}."
            )
        mat = mat.reshape((n, n))

    # Self-regulation is fixed to -1 in the supplied simulation scripts.
    np.fill_diagonal(mat, -1.0)

    return {
        "taxa": np.asarray(taxa, dtype=str),
        "initial": initial.to_numpy(dtype=float).reshape(-1, 1),
        "immigration_weight": immigration.to_numpy(dtype=float).reshape(-1, 1),
        "growth_weight": growth.to_numpy(dtype=float).reshape(-1, 1),
        "interaction": mat,
        "paths": {
            "initial": str(init_path),
            "immigration": str(imm_path),
            "growth": str(gro_path),
            "interaction": str(int_path),
        },
    }


def calibrated_parameters(inputs: Mapping, cfg: Mapping, condition: str,
                          immigration: float, interaction: float,
                          growth: float, noise: float) -> dict:
    noise_multiplier = cfg["warmed_noise_multiplier"] if condition == "warming" else 1.0
    return {
        "initial_condition": np.asarray(inputs["initial"], dtype=float).copy(),
        "immigration_rate": np.asarray(inputs["immigration_weight"], dtype=float) * float(immigration),
        "interaction_matrix": np.asarray(inputs["interaction"], dtype=float).copy(),
        "growth_rate": np.asarray(inputs["growth_weight"], dtype=float) * float(growth),
        "noise_linear": float(noise) * noise_multiplier,
    }


def process_intensity_parameters(inputs: Mapping, cfg: Mapping, condition: str,
                                 immigration: float, interaction: float,
                                 growth: float, noise: float,
                                 process: str, intensity: float) -> dict:
    """
    Build the process-intensity parameterization used in the revised analysis.

    z = 1 is the calibrated model. At z = 0:
      immigration -> 0
      noise       -> 0
      interaction -> off-diagonal interactions 0; diagonal self-regulation -1
      growth      -> homogeneous growth (all taxa receive the scalar g)
    """
    if process not in PROCESSES:
        raise ValueError(f"Unknown process: {process}")
    z = float(intensity)
    if z < 0 or z > 1:
        raise ValueError("Intensity must be in [0, 1].")

    p = calibrated_parameters(inputs, cfg, condition, immigration, interaction, growth, noise)

    if process == "immigration":
        p["immigration_rate"] *= z

    elif process == "noise":
        p["noise_linear"] *= z

    elif process == "interaction":
        full = p["interaction_matrix"]
        scaled = full.copy()
        offdiag = ~np.eye(full.shape[0], dtype=bool)
        scaled[offdiag] = z * full[offdiag]
        np.fill_diagonal(scaled, -1.0)
        p["interaction_matrix"] = scaled

    elif process == "growth":
        species_specific = p["growth_rate"]
        homogeneous = np.ones_like(species_specific) * float(growth)
        p["growth_rate"] = (1.0 - z) * homogeneous + z * species_specific

    return p


def simulate_sglv(params: Mapping, taxa: Sequence[str], seed: int,
                  clip_max: float, dt: float = 1.0,
                  steps_per_snapshot: int = 13,
                  n_snapshots: int = 12) -> pd.DataFrame:
    """
    Simulate the GLV + multiplicative linear Langevin model.

    The Euler--Maruyama update is the one used in the supplied Python scripts:
      dx = (A x * x + immigration + growth*x) dt
           + noise*x*sqrt(dt)*N(0,1)

    The supplied scripts call ``initial=12`` and then use ``range(initial+1)``
    between stored outputs, i.e. 13 integration steps.  This refactor preserves
    that update cadence explicitly while emitting 12 labeled snapshots
    (Y2010--Y2021), which is the shape expected by the supplied resampling code.
    Downstream TPL fitting uses the first 11 snapshots (2010--2020).
    """
    rng = np.random.RandomState(int(seed))
    x = np.asarray(params["initial_condition"], dtype=float).reshape(-1, 1).copy()
    A = np.asarray(params["interaction_matrix"], dtype=float)
    imm = np.asarray(params["immigration_rate"], dtype=float).reshape(-1, 1)
    growth = np.asarray(params["growth_rate"], dtype=float).reshape(-1, 1)
    noise = float(params["noise_linear"])

    n = len(x)
    if A.shape != (n, n):
        raise ValueError(f"Interaction matrix shape {A.shape} does not match N={n}.")

    snapshots = []
    for k in range(int(n_snapshots)):
        for _ in range(int(steps_per_snapshot)):
            dx_det = (A.dot(x) * x + imm + growth * x) * dt
            dx_stoch = noise * x * np.sqrt(dt) * rng.normal(0.0, 1.0, size=x.shape)
            x = x + dx_det + dx_stoch
            x = np.clip(x, 0.0, float(clip_max))
        snapshots.append(x.ravel().copy())

    years = [f"Y{2010+i}" for i in range(int(n_snapshots))]
    out = pd.DataFrame(np.asarray(snapshots), index=years, columns=list(taxa))
    out.index.name = "Year"
    return out


def write_simulation(out: pd.DataFrame, path: Path, overwrite: bool = False) -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        return False
    out.to_csv(path)
    return True


def calibration_grid(cfg: Mapping) -> List[Tuple[float, float, float, int, float]]:
    # Keep interaction strength outermost so only one (potentially very large)
    # interaction matrix needs to reside in memory at a time.
    return [
        (d, a, g, seed, noise)
        for a in cfg["interaction_grid"]
        for d in cfg["immigration_grid"]
        for g in cfg["growth_grid"]
        for seed in SEEDS
        for noise in cfg["noise_grid"]
    ]


def validate_group_names(groups: Sequence[str], configs: Mapping[str, dict]) -> List[str]:
    if len(groups) == 1 and groups[0].lower() == "all":
        return list(configs.keys())
    bad = [g for g in groups if g not in configs]
    if bad:
        raise ValueError(f"Unknown group(s): {bad}. Available: {list(configs)}")
    return list(groups)
