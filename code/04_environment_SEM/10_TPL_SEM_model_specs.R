# ==============================================================================
# Shared model specifications for taxonomic TPL (tTTPL) and phylogenetic TPL
# (pTTPL) mixed-effects piecewise SEMs.
#
# The environmental backbone is intentionally identical between tTTPL and
# pTTPL. Only the taxon-specific beta equations differ between the analyses.
#
# No path selection is performed inside splitplot_piecewiseSEM_engine.R.
# ==============================================================================

make_shared_TPL_environmental_equations <- function() {
  list(
    temperature_mean = c("warm", "drying", "clip"),
    moisture_mean = c("warm", "drying", "wetting", "clip"),
    moisture_cv = c("wetting", "temperature_mean"),
    pH_mean = c("moisture_mean", "moisture_cv"),
    NO3N_cv = c("clip"),
    TC_cv = c("warm", "drying", "moisture_cv"),
    plant_biomass_mean = c("drying", "wetting", "clip"),
    plant_richness_mean = c("warm", "clip", "moisture_mean")
  )
}

make_group_spec <- function(environmental_equations,
                            beta_equations,
                            target,
                            correlated_errors = NULL) {
  out <- list(
    target = target,
    equations = c(
      environmental_equations,
      stats::setNames(list(beta_equations[[target]]), target)
    )
  )
  if (!is.null(correlated_errors) && length(correlated_errors)) {
    out$correlated_errors <- correlated_errors
  }
  out
}


# ==============================================================================
# tTTPL
# ==============================================================================
make_tTTPL_model_spec <- function() {

  environmental_equations <- make_shared_TPL_environmental_equations()

  beta_equations <- list(
    Bac_beta = c(
      "temperature_mean",
      "pH_mean",
      "moisture_cv",
      "TC_cv"
    ),
    Arch_beta = c(
      "moisture_cv",
      "TC_mean",
      "TC_cv"
    ),
    Fun_18S_beta = c(
      "moisture_mean",
      "moisture_cv",
      "NO3N_cv"
    ),
    Pro_beta = c(
      "moisture_mean",
      "plant_richness_mean",
      "moisture_cv",
      "NO3N_cv"
    ),
    Nem_beta = c(
      "moisture_mean",
      "plant_richness_mean",
      "moisture_cv"
    )
  )

  list(
    Bacteria = make_group_spec(
      environmental_equations, beta_equations, "Bac_beta"
    ),
    Archaea = make_group_spec(
      environmental_equations,
      beta_equations,
      "Arch_beta",
      correlated_errors = list(c("TC_mean", "TC_cv"))
    ),
    Fungi = make_group_spec(
      environmental_equations, beta_equations, "Fun_18S_beta"
    ),
    Protists = make_group_spec(
      environmental_equations, beta_equations, "Pro_beta"
    ),
    Nematoda = make_group_spec(
      environmental_equations, beta_equations, "Nem_beta"
    )
  )
}


# ==============================================================================
# pTTPL
# ==============================================================================
make_pTTPL_model_spec <- function() {

  # Keep the full upstream environmental structure identical to tTTPL.
  environmental_equations <- make_shared_TPL_environmental_equations()

  # Current pTTPL endpoint equations. The upstream environmental backbone is
  # identical to tTTPL, while target equations may contain additional exogenous
  # environmental predictors. NO3N_mean and NO3N_cv are retained as distinct
  # variables where explicitly specified; no equivalence between their mean and
  # temporal-variability effects is assumed. No path selection is performed in
  # the piecewise-SEM engine.
  beta_equations <- list(
    Bac_pTTPL_beta = c(
      "temperature_mean",
      "pH_mean"
    ),
    Arch_pTTPL_beta = c(
      "moisture_cv",
      "NO3N_mean"
    ),
    Fun_18S_pTTPL_beta = c(
      "moisture_mean",
      "NO3N_cv"
    ),
    Pro_pTTPL_beta = c(
      "moisture_mean",
      "moisture_cv",
      "plant_richness_mean"
    ),
    Nem_pTTPL_beta = c(
      "moisture_mean",
      "moisture_cv",
      "plant_richness_mean"
    )
  )

  list(
    Bacteria = make_group_spec(
      environmental_equations, beta_equations, "Bac_pTTPL_beta"
    ),
    Archaea = make_group_spec(
      environmental_equations, beta_equations, "Arch_pTTPL_beta",
      correlated_errors = list(c("NO3N_mean", "NO3N_cv"))
    ),
    Fungi = make_group_spec(
      environmental_equations, beta_equations, "Fun_18S_pTTPL_beta"
    ),
    Protists = make_group_spec(
      environmental_equations, beta_equations, "Pro_pTTPL_beta"
    ),
    Nematoda = make_group_spec(
      environmental_equations, beta_equations, "Nem_pTTPL_beta"
    )
  )
}
