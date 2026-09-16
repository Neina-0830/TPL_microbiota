# Data directory

All external inputs used by the canonical analysis workflow are read from `./data/`. The analysis scripts never require a machine-specific absolute path.

## Primary empirical files

| Soil-biotic group | Community table | Taxonomy / trophic annotation | Rooted phylogeny |
|---|---|---|---|
| Archaea | `resample7842_otutab_Archaea_2010to2020.csv` | `resample7842_otutab_Archaea_2010to2020_tax.csv` | `rooted_tree_otus_Archaea_2010to2020_resample7842.nwk` |
| Bacteria | `resample26890_otutab_Bacteria_2010to2020.csv` | `resample26890_otutab_Bacteria_2010to2020_tax.csv` | `rooted_tree_otus_Bacteria_2010to2020_resample26890.nwk` |
| Fungi | `resample541_otutab_Fungi_2010to2020.csv` | `resample541_otutab_Fungi_2010to2020_tax.csv` | `rooted_tree_otus_Fungi_2010to2020_resample541.nwk` |
| Protists | `resample1064_otutab_Protists_2010to2020.csv` | `resample1064_otutab_Protists_2010to2020_tax.csv` | `rooted_tree_otus_Protists_2010to2020_resample1064.nwk` |
| Nematoda | `resample296_otutab_Nematoda_2010to2020.csv` | `resample296_otutab_Nematoda_2010to2020_tax_trophic.csv` | `rooted_tree_otus_Nematoda_2010to2020_resample296.nwk` |

The community tables contain 528 annual samples (48 subplots × 11 years) after group-specific sequencing-depth standardization. Taxa are rows and samples are columns.

Two design/environment files are used across modules:

- `summary of environmental data2010-2020.csv`: sample-level experimental metadata, environmental variables and ecosystem-function measurements;
- `treatment_list10_20.csv`: 48-row subplot-level treatment table.

## Generated prepared files

Run:

```bash
Rscript code/00_prepare_data.R
```

This creates `./data/prepared/`, including harmonized metadata, plot-level relative-abundance matrices, full relative-abundance matrices and taxonomy tables. `./data/prepared/` is generated data and is excluded from version control by `.gitignore`.

## sGLV-specific inputs

The sGLV simulations use additional model-derived parameter files under `./data/sGLV/`. See [`sGLV/README.md`](sGLV/README.md). These files encode the original empirical immigration, growth and interaction structures and are intentionally kept distinct from the primary community tables.
