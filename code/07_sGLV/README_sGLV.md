# Stochastic generalized Lotka–Volterra workflow

Run from the repository root. Parameters are read from `data/sGLV/`; products are written to `results/17_sGLV/`. See `data/sGLV/README.md` for exact filenames and `sGLV_group_config.csv` for group-specific settings.

## 1. Check inputs and calibrate

```bash
python code/07_sGLV/16_check_sGLV_inputs.py
python code/07_sGLV/17_sGLV_calibration_simulate.py --group all --condition both
TPL_SGLV_RESAMPLE_MODE=calibration Rscript code/07_sGLV/18_sGLV_resample.R
Rscript code/07_sGLV/19_sGLV_calibration_TPL.R
Rscript code/07_sGLV/20_sGLV_calibration_community_metrics.R
Rscript code/07_sGLV/21_sGLV_select_top10.R
```

Calibration uses TPL exponents, richness and compositional dissimilarity, with weights specified in `sGLV_calibration_weights.csv`. Candidates are ordered by joint calibration loss. The first ten candidates whose simulated exponents do not differ significantly from empirical estimates in either temperature condition are retained separately for each group and TPL type. The comparison uses two-sided Wilcoxon rank-sum tests. A non-significant difference is a selection criterion, not evidence of model equivalence. Selection stops if fewer than ten eligible candidates exist.

## 2. Compare full and process-removal models

```bash
python code/07_sGLV/22_sGLV_process_removal_simulate.py --tpl tTTPL --group all --condition both
python code/07_sGLV/22_sGLV_process_removal_simulate.py --tpl pTTPL --group all --condition both
TPL_SGLV_RESAMPLE_MODE=process_removal Rscript code/07_sGLV/18_sGLV_resample.R
Rscript code/07_sGLV/23_sGLV_process_removal_TPL.R
Rscript code/07_sGLV/24_sGLV_process_importance.R
```

The five scenarios are `full`, `immigration_removed`, `growth_removed`, `interaction_removed` and `noise_removed`. Full and removal models share each initial condition, calibrated parameter set and stochastic seed. Growth-rate removal homogenizes the calibrated taxon-specific rates to their arithmetic mean. Removing interactions retains self-regulation of −1. Remaining components are held fixed without recalibration.

Three seeds (0, 1, 2) are simulated for each parameter set, condition and scenario. Resampling uses a common seed for matched scenarios; this prevents file-list order from changing the sampling design. Each scenario supplies 30 trajectories per group × TPL type × condition. Temporal slopes use 2010–2020 (the first eleven stored snapshots).

Seed-level beta estimates are averaged within each parameter set before computing effects. Relative importance is the absolute full-minus-removal difference divided by the sum of absolute differences over all four processes. Effects, signed importance, relative importance and parameter-set variation are all retained. There is one full-model baseline per parameter set, condition and seed.

## 3. Main outputs

| File | Contents |
|---|---|
| `Calibration/<TPL>/<Group>_top10_parameters.csv` | Retained parameter sets and calibration diagnostics |
| `Taylor/process_removal_top10/TPL_estimates.csv` | Full and removal exponents for all completed panels |
| `process_importance/process_effects_by_parameter.csv` | Signed effects and normalized absolute/signed importance |
| `process_importance/process_importance_summary.csv` | Ensemble means, bootstrap intervals and effect directions |
| `process_importance/Kruskal_Wallis.csv` | Four-process omnibus comparisons |
| `process_importance/Dunn_Bonferroni.csv` | Six adjusted pairwise comparisons per panel |
| `process_importance/*_process_importance.pdf` | Mean importance, individual parameter sets and bootstrap intervals |

The random seeds are computational realizations. The ten parameter sets characterize model dependence; neither represents independent biological replication. Simulations evaluate conditional model effects and do not establish that immigration increased under experimental warming. Noise scaling uses iCAMP-derived information, so this analysis is not an independent validation of iCAMP.

## Input availability and demo

The model-derived parameter files are not included in the downloaded public repository. Add the original files before running the field-parameterized simulations; the preflight script reports their absence explicitly. Do not substitute synthetic files for these inputs when reproducing the manuscript.

For a self-contained example, run `python code/07_sGLV/demo_process_removal.py`. See the root README for expected outputs and statistical-demo commands.
