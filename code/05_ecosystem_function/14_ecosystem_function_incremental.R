# ==============================================================================
# TPL beta - ecosystem function incremental explanatory-power analysis
# ==============================================================================
# Purpose
# -------
# Directly address the reviewer question:
# "Does the TPL slope add explanatory power over temperature?"
#
# This script analyzes:
#   1) whole-community TPL beta for Bacteria, Archaea, Fungi, Protists, Nematoda
#   2) all phylum / trophic-group beta columns in:
#        All_taxa_diff_Phylum_Trophic_tTTPL_beta_all.csv
#        All_taxa_diff_Phylum_Trophic_pTTPL_beta_all.csv
#   3) both tTTPL and pTTPL
#   4) 12 ecosystem-function properties:
#        GPP_mean, GPP_cv, NEE_mean, NEE_cv, ER_mean, ER_cv,
#        Ra_mean, Ra_cv, Rh_mean, Rh_cv, Rs_mean, Rs_cv
#
# Primary model comparison
# ------------------------
# Base:
#   Function_z ~ Temperature_z + warm + drying + wetting + clip
#                + (1 | block/WholePlotID)
#
# Full:
#   Function_z ~ Temperature_z + Beta_z + warm + drying + wetting + clip
#                + (1 | block/WholePlotID)
#
# The primary quantities are:
#   - standardized conditional coefficient for Beta_z
#   - Kenward-Roger P value for Beta_z
#   - marginal R2 for base and full models
#   - Delta marginal R2 after adding beta
#   - nested-model likelihood-ratio test (ML fits)
#   - Delta AIC
#   - 999 whole-plot cluster-bootstrap CIs for Beta_z and Delta R2m
#
# Sensitivity model comparison
# ----------------------------
# Base:
#   Function_z ~ Temperature_z + (1 | block/WholePlotID)
#
# Full:
#   Function_z ~ Temperature_z + Beta_z + (1 | block/WholePlotID)
#
# This directly matches the reviewer's wording "over temperature".
# By default the sensitivity models are fit without bootstrap because the
# primary analysis already uses 999 whole-plot bootstrap replicates.
#
# IMPORTANT interpretation
# ------------------------
# This is an incremental explanatory-power / conditional-association analysis,
# not a causal effect analysis and not classical Venn-diagram variance
# partitioning. A positive Delta R2m indicates additional fixed-effect variance
# explained after beta is added to the specified base model.
# ==============================================================================


# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

# ------------------------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------------------------
# Run this script from the repository root, or set TPL_PROJECT_ROOT to the
# project directory that contains the data folders documented in README.md.
PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ENV_FILE <-
  "./results/10_environment_SEM_inputs/all_environment_TPL_function_48plots.csv"

TTTPL_SUBGROUP_FILE <-
  "./results/10_environment_SEM_inputs/tTTPL_subgroup_beta_wide.csv"

PTTPL_SUBGROUP_FILE <-
  "./results/10_environment_SEM_inputs/pTTPL_subgroup_beta_wide.csv"

# Primary treatment file used in your existing scripts.
INFO_FILE <-
  "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"

OUTPUT_DIR <-
  "./results/13_ecosystem_function_incremental"

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

# Reviewer-facing ecosystem-function variables
FUNCTION_LIST <- c(
  "GPP_mean", "GPP_cv",
  "NEE_mean", "NEE_cv",
  "ER_mean",  "ER_cv",
  "Ra_mean",  "Ra_cv",
  "Rh_mean",  "Rh_cv",
  "Rs_mean",  "Rs_cv"
)

TEMPERATURE_VAR <- "temperature_mean"

TREATMENT_VARS <- c(
  "Warm",
  "Precipitation",
  "Clip"
)

# Final whole-community variables used in the revised manuscript
WHOLE_TTTPL_TARGETS <- c(
  "Bac_beta",
  "Arch_beta",
  "Fun_18S_beta",
  "Pro_beta",
  "Nem_beta"
)

WHOLE_PTTPL_TARGETS <- c(
  "Bac_pTTPL_beta",
  "Arch_pTTPL_beta",
  "Fun_18S_pTTPL_beta",
  "Pro_pTTPL_beta",
  "Nem_pTTPL_beta"
)

# The older files/scripts may also contain Fun_beta / Fun_pTTPL_beta.
# Keep FALSE for the final manuscript if Fun_18S is the fungal dataset used.
INCLUDE_LEGACY_FUNGAL_COLUMNS <- FALSE

if (INCLUDE_LEGACY_FUNGAL_COLUMNS) {
  WHOLE_TTTPL_TARGETS <-
    c(
      WHOLE_TTTPL_TARGETS,
      "Fun_beta"
    )

  WHOLE_PTTPL_TARGETS <-
    c(
      WHOLE_PTTPL_TARGETS,
      "Fun_pTTPL_beta"
    )
}

# Bootstrap settings
RUN_PRIMARY_BOOTSTRAP <- TRUE
N_BOOT <- 999

# For a quick code check before the final run, temporarily use:
#N_BOOT <- 20

# Parallel bootstrap is supported on macOS/Linux.
# Set N_CORES <- 1 if you prefer serial execution.
N_CORES <- max(
  1,
  min(
    6,
    parallel::detectCores() - 1
  )
)

# Temperature-only sensitivity analysis
RUN_TEMPERATURE_ONLY_SENSITIVITY <- TRUE

# Bootstrap the sensitivity analysis too?
# Usually not necessary for the reviewer response because the primary analysis
# already has whole-plot bootstrap uncertainty.
BOOTSTRAP_TEMPERATURE_ONLY <- FALSE

set.seed(20260902)


# ==============================================================================
# 1. PACKAGES
# ==============================================================================

required_packages <- c(
  "lme4",
  "lmerTest",
  "MuMIn",
  "dplyr",
  "ieggr"
)

missing_packages <-
  required_packages[
    !sapply(
      required_packages,
      requireNamespace,
      quietly = TRUE
    )
  ]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running the script:\n",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)
library(ieggr)


# ==============================================================================
# 2. BASIC HELPER FUNCTIONS
# ==============================================================================

zscore <- function(x) {

  x <- as.numeric(x)

  sx <- sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(sx) ||
    sx == 0
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  as.numeric(
    scale(x)
  )
}


find_info_file <- function() {

  if (
    file.exists(
      INFO_FILE
    )
  ) {
    return(
      INFO_FILE
    )
  }

  if (
    file.exists(
      INFO_FILE
    )
  ) {
    return(
      INFO_FILE
    )
  }

  stop(
    "Could not find the treatment information file.\n",
    "Checked:\n",
    INFO_FILE,
    "\n",
    INFO_FILE
  )
}


align_by_rownames <- function(
    x,
    ref,
    object_name = "object"
) {

  if (
    is.null(
      rownames(x)
    ) ||
    is.null(
      rownames(ref)
    )
  ) {
    stop(
      object_name,
      " or reference data do not have row names."
    )
  }

  missing_rows <-
    setdiff(
      rownames(ref),
      rownames(x)
    )

  if (
    length(
      missing_rows
    ) > 0
  ) {
    stop(
      object_name,
      " is missing ",
      length(missing_rows),
      " reference subplot(s):\n",
      paste(
        missing_rows,
        collapse = ", "
      )
    )
  }

  x[
    rownames(ref),
    ,
    drop = FALSE
  ]
}


safe_r2 <- function(model) {

  ans <-
    try(
      suppressWarnings(
        as.data.frame(
          MuMIn::r.squaredGLMM(
            model
          )
        )
      ),
      silent = TRUE
    )

  if (
    inherits(
      ans,
      "try-error"
    ) ||
    nrow(ans) == 0
  ) {
    return(
      c(
        R2m = NA_real_,
        R2c = NA_real_
      )
    )
  }

  c(
    R2m = as.numeric(
      ans$R2m[1]
    ),
    R2c = as.numeric(
      ans$R2c[1]
    )
  )
}


get_beta_summary <- function(model) {

  sm <-
    coef(
      summary(
        model
      )
    )

  if (
    !"Beta_z" %in%
      rownames(sm)
  ) {
    return(
      c(
        Beta = NA_real_,
        SE = NA_real_,
        P_Satterthwaite = NA_real_
      )
    )
  }

  p_col <-
    grep(
      "^Pr\\(",
      colnames(sm),
      value = TRUE
    )

  p_value <-
    if (
      length(
        p_col
      ) > 0
    ) {
      as.numeric(
        sm[
          "Beta_z",
          p_col[1]
        ]
      )
    } else {
      NA_real_
    }

  c(
    Beta =
      as.numeric(
        sm[
          "Beta_z",
          "Estimate"
        ]
      ),

    SE =
      as.numeric(
        sm[
          "Beta_z",
          "Std. Error"
        ]
      ),

    P_Satterthwaite =
      p_value
  )
}


get_KR_p <- function(model) {

  aa <-
    try(
      suppressMessages(
        lmerTest::anova(
          model,
          type = 2,
          ddf = "Kenward-Roger"
        )
      ),
      silent = TRUE
    )

  if (
    inherits(
      aa,
      "try-error"
    ) ||
    !"Beta_z" %in%
      rownames(aa)
  ) {
    return(
      NA_real_
    )
  }

  p_col <-
    grep(
      "^Pr\\(",
      colnames(aa),
      value = TRUE
    )

  if (
    length(
      p_col
    ) == 0
  ) {
    return(
      NA_real_
    )
  }

  as.numeric(
    aa[
      "Beta_z",
      p_col[1]
    ]
  )
}


bootstrap_two_sided_p <- function(x) {

  x <-
    x[
      is.finite(x)
    ]

  n <-
    length(x)

  if (
    n == 0
  ) {
    return(
      NA_real_
    )
  }

  p_lower <-
    (
      sum(
        x <= 0
      ) + 1
    ) /
    (
      n + 1
    )

  p_upper <-
    (
      sum(
        x >= 0
      ) + 1
    ) /
    (
      n + 1
    )

  min(
    1,
    2 *
      min(
        p_lower,
        p_upper
      )
  )
}


safe_quantile_95 <- function(x) {

  x <-
    x[
      is.finite(x)
    ]

  if (
    length(x) < 20
  ) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }

  qq <-
    quantile(
      x,
      probs = c(
        0.025,
        0.975
      ),
      na.rm = TRUE,
      names = FALSE
    )

  c(
    low =
      as.numeric(
        qq[1]
      ),

    high =
      as.numeric(
        qq[2]
      )
  )
}


# ==============================================================================
# 3. READ AND PREPARE DATA
# ==============================================================================

env <-
  lazyopen(
    ENV_FILE
  )

tttpl_subgroup <-
  lazyopen(
    TTTPL_SUBGROUP_FILE
  )

pttpl_subgroup <-
  lazyopen(
    PTTPL_SUBGROUP_FILE
  )

info_file <-
  find_info_file()

info <-
  lazyopen(
    info_file
  )

# Align every table to the 48-row environmental table
info <-
  align_by_rownames(
    info,
    env,
    "info"
  )

tttpl_subgroup <-
  align_by_rownames(
    tttpl_subgroup,
    env,
    "tTTPL subgroup beta table"
  )

pttpl_subgroup <-
  align_by_rownames(
    pttpl_subgroup,
    env,
    "pTTPL subgroup beta table"
  )


# Check required environmental/function variables
required_env_vars <-
  unique(
    c(
      FUNCTION_LIST,
      TEMPERATURE_VAR,
      WHOLE_TTTPL_TARGETS,
      WHOLE_PTTPL_TARGETS
    )
  )

# Whole-community target variables are checked separately because optional
# legacy fungal columns may not be present.
required_non_target_vars <-
  unique(
    c(
      FUNCTION_LIST,
      TEMPERATURE_VAR
    )
  )

missing_env <-
  setdiff(
    required_non_target_vars,
    colnames(env)
  )

if (
  length(
    missing_env
  ) > 0
) {
  stop(
    "Missing required environmental/function variables in ",
    ENV_FILE,
    ":\n",
    paste(
      missing_env,
      collapse = ", "
    )
  )
}


# Check treatment/design variables
required_info_vars <-
  c(
    "block",
    "WholePlotID",
    TREATMENT_VARS
  )

missing_info <-
  setdiff(
    required_info_vars,
    colnames(info)
  )

if (
  length(
    missing_info
  ) > 0
) {
  stop(
    "Missing treatment/design variables in ",
    info_file,
    ":\n",
    paste(
      missing_info,
      collapse = ", "
    )
  )
}


# Design variables
design_info <-
  info[
    ,
    required_info_vars,
    drop = FALSE
  ]

design_info$block <-
  factor(
    design_info$block
  )

design_info$WholePlotID <-
  factor(
    design_info$WholePlotID
  )
design_info$Warm <- factor(
  design_info$Warm,
  levels = c("unwarming", "warming")
)

design_info$Precipitation <- factor(
  design_info$Precipitation,
  levels = c("normal", "half", "double")
)

design_info$Clip <- factor(
  design_info$Clip,
  levels = c("unclipping", "clipping")
)

# Add subgroup beta columns to separate tTTPL and pTTPL datasets.
# Drop duplicate names if a subgroup file happens to contain a column already
# present in the base environmental table.

tttpl_subgroup_cols <-
  setdiff(
    colnames(
      tttpl_subgroup
    ),
    colnames(env)
  )

pttpl_subgroup_cols <-
  setdiff(
    colnames(
      pttpl_subgroup
    ),
    colnames(env)
  )


dat_tTTPL <-
  cbind(
    env,
    tttpl_subgroup[
      ,
      tttpl_subgroup_cols,
      drop = FALSE
    ],
    design_info
  )


dat_pTTPL <-
  cbind(
    env,
    pttpl_subgroup[
      ,
      pttpl_subgroup_cols,
      drop = FALSE
    ],
    design_info
  )


# ==============================================================================
# 4. DEFINE WHOLE-COMMUNITY AND SUBGROUP TARGETS
# ==============================================================================

whole_t_targets <-
  intersect(
    WHOLE_TTTPL_TARGETS,
    colnames(
      dat_tTTPL
    )
  )

whole_p_targets <-
  intersect(
    WHOLE_PTTPL_TARGETS,
    colnames(
      dat_pTTPL
    )
  )


missing_whole_t <-
  setdiff(
    WHOLE_TTTPL_TARGETS,
    whole_t_targets
  )

missing_whole_p <-
  setdiff(
    WHOLE_PTTPL_TARGETS,
    whole_p_targets
  )

if (
  length(
    missing_whole_t
  ) > 0
) {
  warning(
    "Whole-community tTTPL targets not found and skipped: ",
    paste(
      missing_whole_t,
      collapse = ", "
    )
  )
}

if (
  length(
    missing_whole_p
  ) > 0
) {
  warning(
    "Whole-community pTTPL targets not found and skipped: ",
    paste(
      missing_whole_p,
      collapse = ", "
    )
  )
}


# Every column supplied by the subgroup beta files is treated as a beta target.
# If the subgroup files contain any non-beta metadata columns, remove them here.

sub_t_targets <-
  intersect(
    colnames(
      tttpl_subgroup
    ),
    colnames(
      dat_tTTPL
    )
  )

sub_p_targets <-
  intersect(
    colnames(
      pttpl_subgroup
    ),
    colnames(
      dat_pTTPL
    )
  )


# Remove any accidental duplicated whole-community target name
sub_t_targets <-
  setdiff(
    sub_t_targets,
    whole_t_targets
  )

sub_p_targets <-
  setdiff(
    sub_p_targets,
    whole_p_targets
  )


# Keep only numeric subgroup beta variables.
# This protects against accidental metadata columns.
sub_t_targets <-
  sub_t_targets[
    sapply(
      dat_tTTPL[
        ,
        sub_t_targets,
        drop = FALSE
      ],
      is.numeric
    )
  ]

sub_p_targets <-
  sub_p_targets[
    sapply(
      dat_pTTPL[
        ,
        sub_p_targets,
        drop = FALSE
      ],
      is.numeric
    )
  ]


# Whole-community labels
whole_t_label_map <-
  c(
    Bac_beta = "Bacteria",
    Arch_beta = "Archaea",
    Fun_18S_beta = "Fungi",
    Pro_beta = "Protists",
    Nem_beta = "Nematoda",
    Fun_beta = "Fungi_legacy"
  )

whole_p_label_map <-
  c(
    Bac_pTTPL_beta = "Bacteria",
    Arch_pTTPL_beta = "Archaea",
    Fun_18S_pTTPL_beta = "Fungi",
    Pro_pTTPL_beta = "Protists",
    Nem_pTTPL_beta = "Nematoda",
    Fun_pTTPL_beta = "Fungi_legacy"
  )


make_target_metadata <- function(
    analysis,
    whole_targets,
    subgroup_targets,
    whole_label_map
) {

  whole_meta <-
    data.frame(
      Analysis = analysis,
      Target = whole_targets,
      Target_label =
        unname(
          whole_label_map[
            whole_targets
          ]
        ),
      Target_level =
        "Whole community",
      stringsAsFactors = FALSE
    )


  subgroup_meta <-
    data.frame(
      Analysis = analysis,
      Target = subgroup_targets,
      Target_label = subgroup_targets,
      Target_level =
        "Phylum/Trophic group",
      stringsAsFactors = FALSE
    )


  bind_rows(
    whole_meta,
    subgroup_meta
  )
}


meta_t <-
  make_target_metadata(
    analysis = "tTTPL",
    whole_targets = whole_t_targets,
    subgroup_targets = sub_t_targets,
    whole_label_map = whole_t_label_map
  )

meta_p <-
  make_target_metadata(
    analysis = "pTTPL",
    whole_targets = whole_p_targets,
    subgroup_targets = sub_p_targets,
    whole_label_map = whole_p_label_map
  )


write.csv(
  bind_rows(
    meta_t,
    meta_p
  ),
  file.path(
    OUTPUT_DIR,
    "00_beta_target_metadata.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 5. WHOLE-PLOT CLUSTER BOOTSTRAP
# ==============================================================================

cluster_boot_wholeplot <- function(dat) {

  block_levels <-
    unique(
      as.character(
        dat$block
      )
    )

  boot_blocks <-
    lapply(
      block_levels,
      function(bb) {

        db <-
          dat[
            as.character(
              dat$block
            ) == bb,
            ,
            drop = FALSE
          ]

        wholeplots <-
          unique(
            as.character(
              db$WholePlotID
            )
          )

        sampled_wholeplots <-
          sample(
            wholeplots,
            size =
              length(
                wholeplots
              ),
            replace = TRUE
          )

        pieces <-
          lapply(
            seq_along(
              sampled_wholeplots
            ),
            function(ii) {

              tmp <-
                db[
                  as.character(
                    db$WholePlotID
                  ) ==
                    sampled_wholeplots[ii],
                  ,
                  drop = FALSE
                ]

              # New bootstrap whole-plot identity is essential when a cluster
              # is sampled more than once.
              tmp$WholePlotID_boot <-
                paste0(
                  "Block",
                  bb,
                  "_BootWP",
                  ii
                )

              tmp
            }
          )

        bind_rows(
          pieces
        )
      }
    )


  out <-
    bind_rows(
      boot_blocks
    )

  out$block <-
    factor(
      out$block
    )

  out$WholePlotID_boot <-
    factor(
      out$WholePlotID_boot
    )

  out
}


# ==============================================================================
# 6. FIT ONE TARGET x FUNCTION ANALYSIS
# ==============================================================================

fit_one_incremental_analysis <- function(
    dat,
    analysis,
    target,
    target_label,
    target_level,
    function_var,
    adjustment = c(
      "temperature_treatments",
      "temperature_only"
    ),
    B = 0,
    n_cores = 1
) {

  adjustment <-
    match.arg(
      adjustment
    )


  vars_needed <-
    c(
      function_var,
      target,
      TEMPERATURE_VAR,
      "block",
      "WholePlotID"
    )

  if (
    adjustment ==
      "temperature_treatments"
  ) {
    vars_needed <-
      c(
        vars_needed,
        TREATMENT_VARS
      )
  }


  dd <-
    dat[
      ,
      unique(
        vars_needed
      ),
      drop = FALSE
    ]

  dd <-
    dd[
      complete.cases(dd),
      ,
      drop = FALSE
    ]


  # Minimum design check
  if (
    nrow(dd) < 10 ||
    length(
      unique(
        dd$block
      )
    ) < 2 ||
    length(
      unique(
        dd$WholePlotID
      )
    ) < 4
  ) {

    return(
      data.frame(
        Analysis = analysis,
        Target_level = target_level,
        Target = target,
        Target_label = target_label,
        Ecosystem_function = function_var,
        Adjustment = adjustment,
        Fit_status = "Insufficient complete cases/design",
        stringsAsFactors = FALSE
      )
    )
  }


  # Standardize continuous response and predictors within the complete-case set
  dd$Function_z <-zscore(dd[[function_var]])

  dd$Beta_z <- zscore(dd[[target]])

  dd$Temperature_z <- zscore(dd[[TEMPERATURE_VAR]])


  if (
    anyNA(
      dd$Function_z
    ) ||
    anyNA(
      dd$Beta_z
    ) ||
    anyNA(
      dd$Temperature_z
    )
  ) {

    return(
      data.frame(
        Analysis = analysis,
        Target_level = target_level,
        Target = target,
        Target_label = target_label,
        Ecosystem_function = function_var,
        Adjustment = adjustment,
        Fit_status = "Zero variance after complete-case filtering",
        stringsAsFactors = FALSE
      )
    )
  }


  if (
    adjustment ==
      "temperature_treatments"
  ) {

    rhs_base <- paste0(
      "Temperature_z + ",
      paste(TREATMENT_VARS, collapse = " * ")
    )

  } else {

    rhs_base <-
      "Temperature_z"
  }


  form_base <-
    as.formula(
      paste0(
        "Function_z ~ ",
        rhs_base,
        " + (1 | block/WholePlotID)"
      )
    )


  form_full <-
    as.formula(
      paste0(
        "Function_z ~ ",
        rhs_base,
        " + Beta_z",
        " + (1 | block/WholePlotID)"
      )
    )


  control_full <-
    lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl =
        list(
          maxfun = 200000
        ),
      check.conv.singular =
        "ignore"
    )


  # --------------------------------------------------------------------------
  # Full-data fits
  # --------------------------------------------------------------------------

  # ML models for nested-model comparison / R2 comparison
  m0_ML <-
    try(
      lmerTest::lmer(
        form_base,
        data = dd,
        REML = FALSE,
        control = control_full
      ),
      silent = TRUE
    )

  m1_ML <-
    try(
      lmerTest::lmer(
        form_full,
        data = dd,
        REML = FALSE,
        control = control_full
      ),
      silent = TRUE
    )


  # REML full model for coefficient inference
  m1_REML <-
    try(
      lmerTest::lmer(
        form_full,
        data = dd,
        REML = TRUE,
        control = control_full
      ),
      silent = TRUE
    )


  if (
    inherits(
      m0_ML,
      "try-error"
    ) ||
    inherits(
      m1_ML,
      "try-error"
    ) ||
    inherits(
      m1_REML,
      "try-error"
    )
  ) {

    return(
      data.frame(
        Analysis = analysis,
        Target_level = target_level,
        Target = target,
        Target_label = target_label,
        Ecosystem_function = function_var,
        Adjustment = adjustment,
        Fit_status = "Full-data model failed",
        stringsAsFactors = FALSE
      )
    )
  }


  beta_summary <-
    get_beta_summary(
      m1_REML
    )

  beta_KR_P <-
    get_KR_p(
      m1_REML
    )


  R2_base <-
    safe_r2(
      m0_ML
    )

  R2_full <-
    safe_r2(
      m1_ML
    )


  delta_R2m <-
    unname(
      R2_full[
        "R2m"
      ] -
        R2_base[
          "R2m"
        ]
    )


  # Likelihood-ratio test: models already fit by ML
  lrt <-
    try(
      anova(
        m0_ML,
        m1_ML,
        refit = FALSE
      ),
      silent = TRUE
    )


  LRT_P <-
    if (
      inherits(
        lrt,
        "try-error"
      )
    ) {
      NA_real_
    } else {

      p_col <-
        grep(
          "^Pr\\(",
          colnames(lrt),
          value = TRUE
        )

      if (
        length(
          p_col
        ) == 0
      ) {
        NA_real_
      } else {
        as.numeric(
          lrt[
            2,
            p_col[1]
          ]
        )
      }
    }


  AIC_base <-
    AIC(
      m0_ML
    )

  AIC_full <-
    AIC(
      m1_ML
    )

  # Positive Delta_AIC_beta means that adding beta lowers AIC
  delta_AIC_beta <-
    AIC_base -
    AIC_full


  # --------------------------------------------------------------------------
  # Bootstrap
  # --------------------------------------------------------------------------

  boot_beta <-
    rep(
      NA_real_,
      B
    )

  boot_delta_R2m <-
    rep(
      NA_real_,
      B
    )


  if (
    B > 0
  ) {

    one_boot <- function(bb) {

      db <-
        cluster_boot_wholeplot(
          dd
        )


      # Re-standardize continuous quantities within each bootstrap replicate
      db$Function_z <-zscore(db[[function_var]])
      db$Beta_z <- zscore(db[[target]])
      db$Temperature_z <-zscore(db[[TEMPERATURE_VAR]])

      if (
        anyNA(
          db$Function_z
        ) ||
        anyNA(
          db$Beta_z
        ) ||
        anyNA(
          db$Temperature_z
        )
      ) {
        return(
          c(
            beta = NA_real_,
            delta_R2m = NA_real_
          )
        )
      }


      if (
        adjustment ==
          "temperature_treatments"
      ) {

        rhs_boot <- paste0(
          "Temperature_z + ",
          paste(TREATMENT_VARS, collapse = " * ")
        )

      } else {

        rhs_boot <-
          "Temperature_z"
      }


      f0_boot <-
        as.formula(
          paste0(
            "Function_z ~ ",
            rhs_boot,
            " + (1 | block/WholePlotID_boot)"
          )
        )


      f1_boot <-
        as.formula(
          paste0(
            "Function_z ~ ",
            rhs_boot,
            " + Beta_z",
            " + (1 | block/WholePlotID_boot)"
          )
        )


      mb0 <-
        try(
          suppressWarnings(
            lme4::lmer(
              f0_boot,
              data = db,
              REML = FALSE,
              control = control_full
            )
          ),
          silent = TRUE
        )

      mb1 <-
        try(
          suppressWarnings(
            lme4::lmer(
              f1_boot,
              data = db,
              REML = FALSE,
              control = control_full
            )
          ),
          silent = TRUE
        )


      if (
        inherits(
          mb0,
          "try-error"
        ) ||
        inherits(
          mb1,
          "try-error"
        )
      ) {

        return(
          c(
            beta = NA_real_,
            delta_R2m = NA_real_
          )
        )
      }


      beta_b <-
        if (
          "Beta_z" %in%
            names(
              fixef(
                mb1
              )
            )
        ) {
          as.numeric(
            fixef(
              mb1
            )[
              "Beta_z"
            ]
          )
        } else {
          NA_real_
        }


      r20 <-
        safe_r2(
          mb0
        )[
          "R2m"
        ]

      r21 <-
        safe_r2(
          mb1
        )[
          "R2m"
        ]


      c(
        beta =
          beta_b,

        delta_R2m =
          as.numeric(
            r21 -
              r20
          )
      )
    }


    # macOS/Linux: parallel fork bootstrap
    # Windows: serial fallback
    if (
      .Platform$OS.type !=
        "windows" &&
      n_cores > 1
    ) {

      boot_list <-
        parallel::mclapply(
          X =
            seq_len(B),
          FUN =
            one_boot,
          mc.cores =
            n_cores,
          mc.set.seed =
            TRUE
        )

    } else {

      boot_list <-
        lapply(
          seq_len(B),
          one_boot
        )
    }


    boot_mat <-
      do.call(
        rbind,
        boot_list
      )


    boot_beta <-
      as.numeric(
        boot_mat[
          ,
          "beta"
        ]
      )

    boot_delta_R2m <-
      as.numeric(
        boot_mat[
          ,
          "delta_R2m"
        ]
      )
  }


  beta_CI <-
    if (
      B > 0
    ) {
      safe_quantile_95(
        boot_beta
      )
    } else {
      c(
        low = NA_real_,
        high = NA_real_
      )
    }


  delta_R2_CI <-
    if (
      B > 0
    ) {
      safe_quantile_95(
        boot_delta_R2m
      )
    } else {
      c(
        low = NA_real_,
        high = NA_real_
      )
    }


  beta_boot_mean <-
    if (
      B > 0
    ) {
      mean(
        boot_beta,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }


  delta_R2_boot_mean <-
    if (
      B > 0
    ) {
      mean(
        boot_delta_R2m,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }


  if (
    !is.finite(
      beta_boot_mean
    )
  ) {
    beta_boot_mean <-
      NA_real_
  }

  if (
    !is.finite(
      delta_R2_boot_mean
    )
  ) {
    delta_R2_boot_mean <-
      NA_real_
  }


  # Function labels for supplementary tables
  function_process <-
    sub(
      "_(mean|cv)$",
      "",
      function_var
    )

  function_property <-
    if (
      grepl(
        "_mean$",
        function_var
      )
    ) {
      "Mean"
    } else if (
      grepl(
        "_cv$",
        function_var
      )
    ) {
      "CV"
    } else {
      NA_character_
    }


  data.frame(

    Analysis =
      analysis,

    Target_level =
      target_level,

    Target =
      target,

    Target_label =
      target_label,

    Ecosystem_function =
      function_var,

    Function_process =
      function_process,

    Function_property =
      function_property,

    Adjustment =
      adjustment,

    N_subplots =
      nrow(dd),

    N_wholeplots =
      length(
        unique(
          dd$WholePlotID
        )
      ),

    N_blocks =
      length(
        unique(
          dd$block
        )
      ),


    # Standardized conditional association
    Beta_std =
      unname(
        beta_summary[
          "Beta"
        ]
      ),

    Beta_SE =
      unname(
        beta_summary[
          "SE"
        ]
      ),

    Beta_P_Satterthwaite =
      unname(
        beta_summary[
          "P_Satterthwaite"
        ]
      ),

    Beta_P_KR =
      beta_KR_P,


    # Incremental explanatory power
    R2m_base =
      unname(
        R2_base[
          "R2m"
        ]
      ),

    R2m_full =
      unname(
        R2_full[
          "R2m"
        ]
      ),

    Delta_R2m_beta =
      delta_R2m,


    R2c_base =
      unname(
        R2_base[
          "R2c"
        ]
      ),

    R2c_full =
      unname(
        R2_full[
          "R2c"
        ]
      ),


    # Nested model comparison
    AIC_base =
      AIC_base,

    AIC_full =
      AIC_full,

    Delta_AIC_beta =
      delta_AIC_beta,

    LRT_P =
      LRT_P,


    # Bootstrap beta uncertainty
    Beta_boot_mean =
      beta_boot_mean,

    Beta_boot_CI_low =
      unname(
        beta_CI[
          "low"
        ]
      ),

    Beta_boot_CI_high =
      unname(
        beta_CI[
          "high"
        ]
      ),

    Beta_boot_P =
      if (
        B > 0
      ) {
        bootstrap_two_sided_p(
          boot_beta
        )
      } else {
        NA_real_
      },


    # Bootstrap Delta R2 uncertainty
    Delta_R2m_boot_mean =
      delta_R2_boot_mean,

    Delta_R2m_boot_CI_low =
      unname(
        delta_R2_CI[
          "low"
        ]
      ),

    Delta_R2m_boot_CI_high =
      unname(
        delta_R2_CI[
          "high"
        ]
      ),

    Delta_R2m_boot_P =
      if (
        B > 0
      ) {
        bootstrap_two_sided_p(
          boot_delta_R2m
        )
      } else {
        NA_real_
      },


    N_boot_requested =
      B,

    N_boot_valid_beta =
      sum(
        is.finite(
          boot_beta
        )
      ),

    N_boot_valid_Delta_R2m =
      sum(
        is.finite(
          boot_delta_R2m
        )
      ),

    Singular_base =
      lme4::isSingular(
        m0_ML,
        tol = 1e-4
      ),

    Singular_full =
      lme4::isSingular(
        m1_ML,
        tol = 1e-4
      ),

    Fit_status =
      "OK",

    stringsAsFactors =
      FALSE
  )
}


# ==============================================================================
# 7. RUN A COMPLETE ANALYSIS SET
# ==============================================================================

run_analysis_set <- function(
    dat,
    target_meta,
    adjustment,
    B
) {

  result_list <-
    list()

  counter <-
    1


  for (
    tt in
      seq_len(
        nrow(
          target_meta
        )
      )
  ) {

    target <-
      target_meta$Target[
        tt
      ]

    target_label <-
      target_meta$Target_label[
        tt
      ]

    target_level <-
      target_meta$Target_level[
        tt
      ]

    analysis <-
      target_meta$Analysis[
        tt
      ]


    for (
      ff in
        FUNCTION_LIST
    ) {

      cat(
        "\n",
        adjustment,
        " | ",
        analysis,
        " | ",
        target_level,
        " | ",
        target,
        " | ",
        ff,
        "\n",
        sep = ""
      )


      result_list[[counter]] <-
        fit_one_incremental_analysis(
          dat =
            dat,
          analysis =
            analysis,
          target =
            target,
          target_label =
            target_label,
          target_level =
            target_level,
          function_var =
            ff,
          adjustment =
            adjustment,
          B =
            B,
          n_cores =
            N_CORES
        )


      counter <-
        counter + 1
    }
  }


  bind_rows(
    result_list
  )
}


# ==============================================================================
# 8. PRIMARY ANALYSIS:
#    temperature + treatments, with whole-plot bootstrap
# ==============================================================================

B_primary <-
  if (
    RUN_PRIMARY_BOOTSTRAP
  ) {
    N_BOOT
  } else {
    0
  }


primary_t <-
  run_analysis_set(
    dat =
      dat_tTTPL,
    target_meta =
      meta_t,
    adjustment =
      "temperature_treatments",
    B =
      B_primary
  )


primary_p <-
  run_analysis_set(
    dat =
      dat_pTTPL,
    target_meta =
      meta_p,
    adjustment =
      "temperature_treatments",
    B =
      B_primary
  )


primary_all <-
  bind_rows(
    primary_t,
    primary_p
  )


# ==============================================================================
# 9. MULTIPLE-TESTING CORRECTION
# ==============================================================================

# Scientific family:
# the 12 ecosystem-function properties tested for one beta target.
#
# Therefore Holm correction is applied separately for each:
#   Analysis x Target
#
# This means, for example, Bacteria tTTPL beta has one family of 12 tests,
# while a bacterial phylum-level beta has its own family of 12 tests.

primary_all <-
  primary_all %>%
  group_by(
    Analysis,
    Target_level,
    Target
  ) %>%
  mutate(

    Beta_P_KR_Holm =
      p.adjust(
        Beta_P_KR,
        method = "holm"
      ),

    LRT_P_Holm =
      p.adjust(
        LRT_P,
        method = "holm"
      ),

    Beta_boot_P_Holm =
      p.adjust(
        Beta_boot_P,
        method = "holm"
      )

  ) %>%
  ungroup()


# Bootstrap-support flags
primary_all <-
  primary_all %>%
  mutate(

    Beta_boot_CI_excludes_zero =
      is.finite(
        Beta_boot_CI_low
      ) &
      is.finite(
        Beta_boot_CI_high
      ) &
      (
        Beta_boot_CI_low > 0 |
        Beta_boot_CI_high < 0
      ),

    Delta_R2m_boot_CI_above_zero =
      is.finite(
        Delta_R2m_boot_CI_low
      ) &
      Delta_R2m_boot_CI_low > 0,

    Added_information_supported =
      Beta_boot_CI_excludes_zero &
      Delta_R2m_boot_CI_above_zero
  )


write.csv(
  primary_all,
  file.path(
    OUTPUT_DIR,
    "01_PRIMARY_temperature_treatments_plus_beta_all_results.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 10. TEMPERATURE-ONLY SENSITIVITY ANALYSIS
# ==============================================================================

if (
  RUN_TEMPERATURE_ONLY_SENSITIVITY
) {

  B_sensitivity <-
    if (
      BOOTSTRAP_TEMPERATURE_ONLY
    ) {
      N_BOOT
    } else {
      0
    }


  sensitivity_t <-
    run_analysis_set(
      dat =
        dat_tTTPL,
      target_meta =
        meta_t,
      adjustment =
        "temperature_only",
      B =
        B_sensitivity
    )


  sensitivity_p <-
    run_analysis_set(
      dat =
        dat_pTTPL,
      target_meta =
        meta_p,
      adjustment =
        "temperature_only",
      B =
        B_sensitivity
    )


  sensitivity_all <-
    bind_rows(
      sensitivity_t,
      sensitivity_p
    )


  sensitivity_all <-
    sensitivity_all %>%
    group_by(
      Analysis,
      Target_level,
      Target
    ) %>%
    mutate(

      Beta_P_KR_Holm =
        p.adjust(
          Beta_P_KR,
          method = "holm"
        ),

      LRT_P_Holm =
        p.adjust(
          LRT_P,
          method = "holm"
        )

    ) %>%
    ungroup()


  write.csv(
    sensitivity_all,
    file.path(
      OUTPUT_DIR,
      "02_SENSITIVITY_temperature_only_plus_beta_all_results.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

} else {

  sensitivity_all <-
    data.frame()
}


# ==============================================================================
# 11. MANUSCRIPT-READY SUPPLEMENTARY TABLE
# ==============================================================================

supplementary_table <-
  primary_all %>%
  select(

    Analysis,

    Target_level,

    Target,

    Target_label,

    Function_process,

    Function_property,

    Ecosystem_function,

    N_subplots,

    N_wholeplots,

    N_blocks,


    # Conditional beta association
    Beta_std,

    Beta_SE,

    Beta_P_KR,

    Beta_P_KR_Holm,

    Beta_boot_CI_low,

    Beta_boot_CI_high,


    # Incremental fixed-effect variance explained
    R2m_base,

    R2m_full,

    Delta_R2m_beta,

    Delta_R2m_boot_CI_low,

    Delta_R2m_boot_CI_high,


    # Nested-model support
    LRT_P,

    LRT_P_Holm,

    Delta_AIC_beta,


    # Diagnostics/support
    Beta_boot_CI_excludes_zero,

    Delta_R2m_boot_CI_above_zero,

    Added_information_supported,

    Singular_full,

    Fit_status
  )


write.csv(
  supplementary_table,
  file.path(
    OUTPUT_DIR,
    "03_Supplementary_Table_TPL_beta_additional_explanatory_power.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 12. WHOLE-COMMUNITY-ONLY SUPPLEMENTARY TABLE
# ==============================================================================

whole_community_table <-
  supplementary_table %>%
  filter(
    Target_level ==
      "Whole community"
  )


write.csv(
  whole_community_table,
  file.path(
    OUTPUT_DIR,
    "04_Whole_community_only_additional_explanatory_power.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 13. PHYLUM / TROPHIC-GROUP-ONLY SUPPLEMENTARY TABLE
# ==============================================================================

subgroup_table <-
  supplementary_table %>%
  filter(
    Target_level ==
      "Phylum/Trophic group"
  )


write.csv(
  subgroup_table,
  file.path(
    OUTPUT_DIR,
    "05_Phylum_Trophic_group_additional_explanatory_power.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 14. COMPACT SUPPORTED-RESULT TABLE
# ==============================================================================

supported_results <-
  supplementary_table %>%
  filter(
    Fit_status ==
      "OK"
  ) %>%
  filter(
    Beta_P_KR_Holm <= 0.05 |
    LRT_P_Holm <= 0.05 |
    Added_information_supported
  ) %>%
  arrange(
    Analysis,
    Target_level,
    Target,
    Ecosystem_function
  )


write.csv(
  supported_results,
  file.path(
    OUTPUT_DIR,
    "06_Supported_results_for_review_and_text.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 15. OPTIONAL: MERGE PRIMARY AND TEMPERATURE-ONLY SENSITIVITY
# ==============================================================================

if (
  nrow(
    sensitivity_all
  ) > 0
) {

  comparison_table <-
    primary_all %>%
    select(
      Analysis,
      Target_level,
      Target,
      Target_label,
      Ecosystem_function,

      Beta_std_primary =
        Beta_std,

      Beta_P_KR_Holm_primary =
        Beta_P_KR_Holm,

      Delta_R2m_primary =
        Delta_R2m_beta,

      LRT_P_Holm_primary =
        LRT_P_Holm
    ) %>%

    left_join(

      sensitivity_all %>%
        select(
          Analysis,
          Target_level,
          Target,
          Ecosystem_function,

          Beta_std_temperature_only =
            Beta_std,

          Beta_P_KR_Holm_temperature_only =
            Beta_P_KR_Holm,

          Delta_R2m_temperature_only =
            Delta_R2m_beta,

          LRT_P_Holm_temperature_only =
            LRT_P_Holm
        ),

      by =
        c(
          "Analysis",
          "Target_level",
          "Target",
          "Ecosystem_function"
        )
    ) %>%

    mutate(

      Same_beta_direction =
        sign(
          Beta_std_primary
        ) ==
        sign(
          Beta_std_temperature_only
        )
    )


  write.csv(
    comparison_table,
    file.path(
      OUTPUT_DIR,
      "07_Primary_vs_temperature_only_sensitivity.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )
}


# ==============================================================================
# 16. SAVE R OBJECTS
# ==============================================================================

saveRDS(
  primary_all,
  file.path(
    OUTPUT_DIR,
    "primary_all_results.rds"
  )
)

saveRDS(
  sensitivity_all,
  file.path(
    OUTPUT_DIR,
    "temperature_only_sensitivity_results.rds"
  )
)


# ==============================================================================
# 17. CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "Analysis completed.\n"
)

cat(
  "Output directory:\n",
  normalizePath(
    OUTPUT_DIR
  ),
  "\n"
)

cat(
  "Primary result rows: ",
  nrow(
    primary_all
  ),
  "\n",
  sep = ""
)

cat(
  "Whole-community targets: ",
  nrow(
    whole_community_table
  ),
  " target-function rows\n",
  sep = ""
)

cat(
  "Phylum/trophic-group targets: ",
  nrow(
    subgroup_table
  ),
  " target-function rows\n",
  sep = ""
)

cat(
  "Rows flagged for review/text: ",
  nrow(
    supported_results
  ),
  "\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)


print(
  supported_results %>%
    select(
      Analysis,
      Target_level,
      Target_label,
      Target,
      Ecosystem_function,
      Beta_std,
      Beta_P_KR_Holm,
      Delta_R2m_beta,
      Delta_R2m_boot_CI_low,
      Delta_R2m_boot_CI_high,
      LRT_P_Holm,
      Added_information_supported
    )
)


# ==============================================================================
# END
# ==============================================================================
