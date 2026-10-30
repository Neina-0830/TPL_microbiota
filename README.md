# Temporal variance scaling across diverse soil microbiota

Analysis workflow for **Experimental Warming Amplifies Temporal Variance Scaling Across Diverse Soil Biota**.

The experiment comprises four blocks, six whole plots per block and two subplots per whole plot: 24 whole plots and 48 subplot trajectories. Annual soil samples cover 2010–2020 (528 samples). Warming and precipitation were applied to whole plots; clipping was applied to the southern subplot. Five soil biotic groups are analyzed: bacteria, archaea, fungi, protists and nematodes.

- Taxonomic temporal TPL (tTTPL) fits log10 temporal variance against log10 temporal mean relative abundance across taxa within a subplot.
- Phylogenetic temporal TPL (pTTPL) uses normalized taxon-specific contributions to abundance-weighted phylogenetic structure, MPD_i,k.

## Files and directories

| Location | Contents |
|---|---|
| `data/` | Public community tables, taxonomy, rooted trees and experimental metadata |
| `data/sGLV/` | Input specification for the model-derived immigration, growth and interaction structures |
| `code/00_prepare_data.R` | Harmonize metadata and produce subplot-level abundance matrices |
| `code/01_empirical_TPL/` | TPL fitting, detrending, taxon removal and treatment inference |
| `code/02_tipshuffle_null/` | Tip-shuffle nulls and phylogroup scaling/warming comparisons |
| `code/03_subgroup_TPL/` | Phylum and nematode trophic-group analyses |
| `code/04_environment_SEM/` | Environmental summaries and mixed-effects piecewise SEM |
| `code/05_ecosystem_function/` | Standardized associations and variation partitioning |
| `code/06_iCAMP/` | Community assembly-process analysis |
| `code/07_sGLV/` | Calibration and single-process removal simulations |
| `CODE_LIST_GITHUB.csv` | Script inventory |
| `DATA_MANIFEST.csv` | Primary data sizes and SHA256 checksums |
| `RUN_ORDER.txt` | Complete ordered commands |

All commands are run from the repository root. External inputs are read from `data/`; generated files are written to `data/prepared/` and `results/`. `TPL_PROJECT_ROOT` can specify a different root directory for batch execution.

## Software and installation

The manuscript reports R 4.3.2 and Python 3.11.8. This repository does not include a frozen R package environment. Install package releases compatible with your R version; if a current CRAN release requires a newer R version, use a compatible archived release or an appropriately updated environment and record `sessionInfo()`.

```bash
Rscript code/install_R_packages.R
python -m pip install -r code/07_sGLV/requirements_sGLV.txt
Rscript code/validate_R_syntax.R
python code/validate_repository.py
```

Required R packages are enumerated in the installation script. Python simulation code uses NumPy and pandas. Standard CPU hardware is sufficient for the synthetic demo. Full empirical phylogenetic matrices, bootstraps, iCAMP and calibration grids require substantially more memory and time; the bacterial distance matrix is particularly large. Installation time depends on whether compiled R dependencies are already available and has not been benchmarked.

## Small demonstration

```bash
python code/07_sGLV/demo_process_removal.py
```

This creates 300 synthetic tTTPL estimates (10 parameter sets × 2 conditions × 3 seeds × 5 scenarios), and 20 normalized process-importance panels under `results/demo_process_removal/`. It does not use field data or reproduce empirical estimates. The observed runtime in the validation environment was approximately one second; runtime varies by machine.

Expected outputs:

- `TPL_estimates.csv`: one estimate per parameter set, condition, seed and scenario;
- `expected_importance.csv`: four relative effects per parameter set and condition, summing to one.

After installing R dependencies, the complete statistical summary can also be demonstrated with:

```bash
TPL_SGLV_BETA_FILE=results/demo_process_removal/TPL_estimates.csv \
TPL_SGLV_IMPORTANCE_OUT=results/demo_process_removal/statistics \
Rscript code/07_sGLV/24_sGLV_process_importance.R
```

## Empirical workflow

Run the empirical commands in `RUN_ORDER.txt` in order. Data preparation verifies annual sample alignment and writes harmonized treatment metadata. Fungal source filenames use `Fungi`; sGLV uses the group label `Fungi_18S` for the same fungal dataset.

### TPL estimates and treatment effects

Both raw and linearly detrended exponents are reported with slope standard errors. Detrended temporal variances are calculated from residuals while retaining original temporal means. Variance is also evaluated at a common reference mean within each group and formulation.

The primary model is:

```r
metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)
```

Precipitation is a three-level categorical factor. Overall warming, reduced precipitation, increased precipitation and clipping effects are estimated by equal-weight marginal contrasts, with Holm adjustment across four effects within each group, formulation and metric. Eleven treatment-versus-control comparisons form a separate Holm family. Cohen's d is descriptive. Multilevel meta-analytic models incorporate SE(beta)^2 to assess uncertainty in the overall treatment conclusions.

Taxon-removal robustness analyses apply to tTTPL: strongly trending taxa are removed, individual leading taxa are omitted, and random 1%, 5% and 10% subsets are removed in 200 realizations per level. Remaining relative abundances are renormalized. These analyses do not estimate aggregate community stability.

Primary treatment outputs include `overall.csv`, `specific.csv`, `uncertainty.csv`, `overall_effect_sizes_and_inference.csv` and treatment-effect PDFs in the corresponding factorial-statistics directories.

### Phylogenetic analyses

The tip-shuffle null reassigns taxon labels among fixed tree tips while retaining complete abundance trajectories, occurrence patterns, topology and branch lengths. The default is 100 shuffles, followed by paired observed-versus-null Wilcoxon tests across 48 subplots and Holm adjustment across groups. A split-plot model assesses warming-induced changes in observed beta_p relative to the subplot null expectation.

The phylogroup module uses taxa with at least 10% prevalence and a fixed group-specific distance quantile: 0.01 for bacteria, 0.40 for archaea, 0.20 for fungi, 0.10 for protists and 0.40 for nematodes. Bacterial and protistan focal taxa are sampled (1,000 and 500, respectively); all eligible focal taxa are used for the other groups. Group beta estimates must be finite in all 48 subplots. Label permutations test the association of phylogenetic distance with differences in mean beta or block-balanced warming response. P values are Holm-adjusted across five groups separately for the two analyses.

### Environmental pathways

Each mixed-effects piecewise SEM component includes block and nested whole-plot random intercepts. The environmental backbone and group-specific endpoint equations are defined in `10_TPL_SEM_model_specs.R`. Directed path coefficients are standardized; indirect effects are products along directed paths, and total effects combine direct and indirect effects. Residual correlated errors are excluded from those sums. Uncertainty uses 999 within-block whole-plot bootstrap samples; four leave-one-block-out refits assess block dependence. Paths are conditional associations within the specified model.

### Ecosystem-function analyses

`14_environment_function_associations.R` fits standardized mixed-effects associations:

```r
z(Function) ~ z(beta) + (1 | Block/WholePlotID)
```

Function is the response and beta is the predictor. Environmental associations use the corresponding standardized beta endpoint and environmental predictor. Tables report standardized coefficients, Kenward–Roger tests, confidence intervals, raw P values and explicitly labeled Holm adjustments.

`15_temperature_TPL_variation_partitioning.R` partitions each function between mean soil temperature and one TPL exponent at a time, conditioning on block. It retains signed adjusted R² fractions for unique temperature, unique beta, shared and unexplained variation. Outputs distinguish fractions of total variation from fractions of the variation remaining after block adjustment. Negative adjusted fractions are retained.

Permutation tests shuffle whole plots within blocks, preserve the paired clipping positions and use 9,999 permutations. The script reports nominal P values and explicitly defined BH-adjusted families. Complementary two-predictor hierarchical partitioning averages contributions over the two predictor-entry orders. The mean and CV function profiles are also evaluated jointly. Outputs include `variation_fractions.csv`, `permutation_tests.csv`, `hierarchical_partitioning.csv` and `B16_single_beta_summary.csv` in a timestamped run directory.

### Community assembly

The iCAMP script records its settings and writes process summaries under `results/14_iCAMP/`. Its default validated settings are for bacteria. Group-specific settings must be supplied for other groups; bacterial binning settings should not be assumed universal.

## sGLV calibration and process removal

See `code/07_sGLV/README_sGLV.md` and `data/sGLV/README.md`.

The process comparison evaluates the full calibrated model against four scenarios, holding all other parameters fixed without recalibration:

| Scenario | Modification |
|---|---|
| Immigration removed | Set immigration to zero |
| Growth-rate differences removed | Replace taxon-specific growth rates by their mean in that calibrated community |
| Interspecific interactions removed | Set off-diagonal coefficients to zero and retain diagonal coefficients of −1 |
| Stochastic noise removed | Set the stochastic-noise amplitude to zero |

For parameter set j, process k and condition c:

```text
effect[j,k,c] = mean_seed(beta_full[j,c]) − mean_seed(beta_removed[j,k,c])
I[j,k,c]      = abs(effect[j,k,c]) / sum_k(abs(effect[j,k,c]))
```

Three stochastic realizations are averaged before normalizing effects within each of the ten retained parameter sets. The output also retains signed effects and signed normalized effects. All-zero denominators are undefined and are not silently assigned zero importance. Summaries and bootstrap intervals describe variation across parameter sets, not biological replication.

Between-process statistics follow the supplied analysis: Kruskal–Wallis tests followed by Dunn tests with Bonferroni correction across six process pairs within each panel. Compact letters are gated by the omnibus test. Individual parameter-set values are shown in plots. These tests describe the calibrated ensemble; the model does not demonstrate a field causal mechanism or that warming increased dispersal.

**Input availability:** The included community tables support the empirical analyses. The public repository snapshot does not contain the model-derived immigration, growth and interaction input files needed for the complete sGLV calibration. These are not fabricated or inferred by the workflow. The preflight script lists missing files. Full simulation reproduction requires those files and cannot proceed until they are supplied. The synthetic demo remains runnable without them.

## Validation and limits

`VALIDATION.md` records checks performed on this release. Syntax and small synthetic examples are distinct from regenerating all manuscript results. The full R workflow, phylogenetic distance calculations and field-parameterized sGLV ensemble have not been executed in the packaging environment. The parameter inputs and original software environment must be supplied for complete numerical reproduction.

## License and citation

Code is distributed under the MIT License. Cite the accompanying manuscript when using this workflow or the associated data. Public data identifiers and manuscript publication details should accompany a formal archival release.
