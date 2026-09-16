# sGLV model-derived inputs

The main empirical TPL workflow uses the raw community, taxonomy, tree, treatment and environmental files stored directly in `./data/`.

The sGLV simulations additionally require the **model-derived parameter inputs used in the original calibration**. These files should be placed under `./data/sGLV/`; they are deliberately not regenerated from the raw community tables because reconstructing them would require re-running the original neutral-model, copy-number/growth and co-occurrence-network parameterization steps and could alter the fitted model.

For each group (`Archaea`, `Bacteria`, `Fungi_18S`, `Protists`, `Nematoda`) and condition (`warming`, `unwarming`), the expected files are:

```text
data/sGLV/
├── rawdata/
│   ├── resample<depth>_otutab_<Group>_<condition>_2010_mean.csv
│   ├── otutab_<Group>_2010to2020_<condition>_relative_immigration.csv
│   └── otutab_<Group>_2010to2020_<condition>_copynumber_relative_growth.csv
└── interaction_matrix/
    └── data1_data2/
        └── otutab_<Group>_<condition>_mean0_sd<strength>_interaction_matrix.csv
```

The required interaction strengths and group-specific grids are recorded in `code/07_sGLV/sGLV_group_config.csv`.

Run the preflight check before a calibration job:

```bash
python code/07_sGLV/16_check_sGLV_inputs.py
```

All generated simulations, resamples, calibration tables and process-intensity outputs are written under `./results/17_sGLV/`.
