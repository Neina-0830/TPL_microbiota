# ==============================================================================
# Sensitivity engine: multilevel path model using CRAN package plssem
#
# This file requires splitplot_piecewiseSEM_engine.R to be sourced first because
# it reuses data/design QC, DAG validation, path enumeration, and bootstrap
# summary helpers. No piecewiseSEM model is fitted by this engine.
#
# IMPORTANT:
# - The manuscript's primary hierarchy is (1 | block/WholePlotID).
# - The current plssem multilevel estimator accepts only ONE clustering
#   variable in lme4-style syntax. It therefore cannot fit both random
#   intercepts in (1 | block/WholePlotID).
# - This optional sensitivity analysis uses the closest supported adjustment:
#   block is represented by K-1 fixed contrasts and BlockWholePlotID is the
#   single random-intercept cluster.
# - The package's ordinary row bootstrap is NOT used. Optional uncertainty is
#   obtained by resampling complete whole plots within each block.
# - With observed single variables rather than multi-indicator latent constructs,
#   this is a multilevel PLS path-model sensitivity analysis, not a latent
#   measurement model and not a solution to low power at n = 48.
# ==============================================================================

resolve_column <- function(x, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) return(hit[1])
  lower <- tolower(names(x))
  idx <- match(tolower(candidates), lower, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) return(names(x)[idx[1]])
  if (required) {
    stop(
      "Could not find any of these columns in parameter_estimates(): ",
      paste(candidates, collapse = ", "),
      ". Available columns: ", paste(names(x), collapse = ", ")
    )
  }
  NULL
}

add_plssem_supported_grouping <- function(dat) {
  dat$block <- droplevels(factor(dat$block))
  dat$WholePlotID <- droplevels(factor(dat$WholePlotID))
  dat$BlockWholePlotID <- interaction(
    dat$block, dat$WholePlotID,
    drop = TRUE, lex.order = TRUE
  )

  # plssem permits one clustering variable only. Add K-1 numeric block
  # contrasts as nuisance fixed effects so that block heterogeneity is still
  # controlled without asking plslmer() to parse a second random cluster.
  block_mm <- stats::model.matrix(
    ~ block,
    data = dat,
    contrasts.arg = list(block = "contr.treatment")
  )
  block_mm <- block_mm[, colnames(block_mm) != "(Intercept)", drop = FALSE]
  block_covariates <- paste0("block_fixed_", seq_len(ncol(block_mm)))
  if (ncol(block_mm)) {
    colnames(block_mm) <- block_covariates
    for (v in block_covariates) dat[[v]] <- as.numeric(block_mm[, v])
  }

  list(data = dat, block_covariates = block_covariates)
}

make_plssem_syntax <- function(equations, block_covariates) {
  lines <- vapply(names(equations), function(response) {
    rhs <- c(
      unique(equations[[response]]),
      block_covariates,
      "(1 | BlockWholePlotID)"
    )
    paste(response, "~", paste(rhs, collapse = " + "))
  }, character(1))
  paste(lines, collapse = "\n")
}

fit_plssem_group <- function(dat, equations) {
  dd <- add_plssem_supported_grouping(dat)
  syntax <- make_plssem_syntax(equations, dd$block_covariates)
  fit <- plssem::pls(
    syntax = syntax,
    data = dd$data,
    standardize = TRUE,
    consistent = FALSE,
    bootstrap = FALSE,
    missing = "listwise",
    verbose = FALSE
  )
  list(
    fit = fit,
    data = dd$data,
    syntax = syntax,
    block_covariates = dd$block_covariates,
    random_structure = "block fixed contrasts + (1|BlockWholePlotID)"
  )
}

extract_plssem_edges <- function(fit_object, equations) {
  tab <- as.data.frame(plssem::parameter_estimates(fit_object$fit))
  lhs_col <- resolve_column(tab, c("lhs", "LHS", "response"))
  op_col <- resolve_column(tab, c("op", "operator"), required = FALSE)
  rhs_col <- resolve_column(tab, c("rhs", "RHS", "predictor"))
  est_col <- resolve_column(tab, c("est", "estimate", "Estimate"))
  se_col <- resolve_column(tab, c("se", "std.error", "Std.Error", "SE"), required = FALSE)
  p_col <- resolve_column(tab, c("pvalue", "p.value", "P.Value", "p"), required = FALSE)

  if (!is.null(op_col)) tab <- tab[as.character(tab[[op_col]]) == "~", , drop = FALSE]
  rows <- list()
  k <- 1L
  for (response in names(equations)) {
    for (parent in equations[[response]]) {
      ii <- as.character(tab[[lhs_col]]) == response & as.character(tab[[rhs_col]]) == parent
      if (sum(ii) != 1L) {
        stop(
          "Expected one structural estimate for ", parent, " -> ", response,
          "; found ", sum(ii), "."
        )
      }
      estimate <- as.numeric(tab[[est_col]][ii])
      se <- if (is.null(se_col)) NA_real_ else as.numeric(tab[[se_col]][ii])
      p <- if (is.null(p_col)) NA_real_ else as.numeric(tab[[p_col]][ii])
      if (!is.finite(p) && is.finite(se) && se > 0) {
        p <- 2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE)
      }
      rows[[k]] <- data.frame(
        From = parent,
        To = response,
        Std_Estimate = estimate,
        Std_SE = se,
        P_value = p,
        Direction = ifelse(estimate >= 0, "positive", "negative"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

bootstrap_plssem_group <- function(dat, equations, sources, target, n_boot, seed, n_cores) {
  if (n_boot <= 0L) return(NULL)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_boot)

  one <- function(i) {
    set.seed(seeds[i])
    dboot <- resample_wholeplots_within_blocks(dat)
    tryCatch({
      fit <- fit_plssem_group(dboot, equations)
      edges <- extract_plssem_edges(fit, equations)
      paths <- enumerate_path_effects(edges, sources, target)
      totals <- total_effects_from_paths(paths, sources = sources, target = target)
      edges$Replicate <- i
      paths$Replicate <- i
      totals$Replicate <- i
      list(edges = edges, paths = paths, totals = totals, error = NULL)
    }, error = function(e) {
      list(edges = NULL, paths = NULL, totals = NULL, error = conditionMessage(e))
    })
  }

  if (.Platform$OS.type == "windows" || n_cores <= 1L) {
    res <- lapply(seq_len(n_boot), one)
  } else {
    res <- parallel::mclapply(seq_len(n_boot), one, mc.cores = n_cores)
  }
  errors <- vapply(res, function(x) x$error %||% "", character(1))
  valid <- errors == ""
  edges <- if (any(valid)) do.call(rbind, lapply(res[valid], `[[`, "edges")) else NULL
  paths <- if (any(valid)) do.call(rbind, lapply(res[valid], `[[`, "paths")) else NULL
  totals <- if (any(valid)) do.call(rbind, lapply(res[valid], `[[`, "totals")) else NULL

  list(
    edges = edges,
    paths = paths,
    totals = totals,
    edge_summary = summarize_bootstrap(edges, c("From", "To"), "Std_Estimate"),
    path_summary = summarize_bootstrap(paths, c("Source_variable", "Target", "Path"), "Path_effect"),
    direct_summary = summarize_bootstrap(totals, c("Source_variable", "Target"), "Direct_effect"),
    indirect_summary = summarize_bootstrap(totals, c("Source_variable", "Target"), "Indirect_effect"),
    total_summary = summarize_bootstrap(totals, c("Source_variable", "Target"), "Total_effect"),
    failures = data.frame(
      Replicate = which(!valid), Error = errors[!valid],
      stringsAsFactors = FALSE
    )
  )
}

lobo_plssem_group <- function(dat, equations, full_edges) {
  blocks <- levels(dat$block)
  rows <- list()
  k <- 1L
  for (b in blocks) {
    dtrain <- droplevels(dat[dat$block != b, , drop = FALSE])
    ans <- tryCatch({
      fit <- fit_plssem_group(dtrain, equations)
      extract_plssem_edges(fit, equations)
    }, error = function(e) NULL)
    if (is.null(ans)) next
    ans$Held_out_block <- b
    rows[[k]] <- ans
    k <- k + 1L
  }
  if (!length(rows)) return(list(estimates = data.frame(), stability = data.frame()))
  est <- do.call(rbind, rows)
  keys <- unique(est[, c("From", "To"), drop = FALSE])
  stability <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    dd <- est[est$From == keys$From[i] & est$To == keys$To[i], , drop = FALSE]
    full <- full_edges$Std_Estimate[
      full_edges$From == keys$From[i] & full_edges$To == keys$To[i]
    ][1]
    data.frame(
      From = keys$From[i], To = keys$To[i],
      Full_Std_Estimate = full,
      LOBO_min = min(dd$Std_Estimate, na.rm = TRUE),
      LOBO_max = max(dd$Std_Estimate, na.rm = TRUE),
      LOBO_mean = mean(dd$Std_Estimate, na.rm = TRUE),
      Sign_consistency_vs_full = mean(sign(dd$Std_Estimate) == sign(full), na.rm = TRUE),
      N_LOBO_folds = nrow(dd),
      stringsAsFactors = FALSE
    )
  }))
  list(estimates = est, stability = stability)
}

run_splitplot_plssem <- function(
    data_file,
    meta_file,
    out_dir,
    analysis_label,
    figure_label,
    model_spec,
    target_labels,
    model_spec_basis = "fixed paths supplied before plssem fitting",
    treatments = c("warm", "drying", "wetting", "clip"),
    alpha_display = 0.05,
    run_bootstrap = TRUE,
    n_boot = 999L,
    seed = 8302026L,
    n_cores = 1L,
    run_lobo = TRUE) {

  assert_packages("plssem")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "groups"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "combined"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "bootstrap"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "LOBO_robustness"), recursive = TRUE, showWarnings = FALSE)

  input <- read_and_match_splitplot(data_file, meta_file, treatments)
  dat <- input$data
  validate_model_spec(model_spec, names(dat), treatments)

  fixed_paths <- do.call(rbind, lapply(names(model_spec), function(group) {
    x <- spec_edges(model_spec[[group]]$equations)
    x$Group <- group
    x$Target <- model_spec[[group]]$target
    x
  }))
  unique_fixed_edges <- unique(fixed_paths[, c("From", "To"), drop = FALSE])

  write_csv(data.frame(
    Analysis = analysis_label,
    N_subplots = nrow(dat),
    N_blocks = length(unique(dat$block)),
    N_wholeplots = length(unique(interaction(
      dat$block, dat$WholePlotID,
      drop = TRUE, lex.order = TRUE
    ))),
    Random_structure_requested = "(1|block/WholePlotID)",
    Plssem_cluster_limit = "one clustering variable",
    Block_adjustment_fitted = "K-1 fixed block contrasts",
    Random_structure_fitted = "(1|BlockWholePlotID)",
    Model_specification_basis = model_spec_basis,
    N_unique_fixed_edges = nrow(unique_fixed_edges),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "00_design_QC.csv"))
  write_csv(input$pair_qc, file.path(out_dir, "00_wholeplot_pair_QC.csv"))
  write_csv(
    fixed_paths,
    file.path(out_dir, "00_fixed_paths_by_group_before_fitting.csv")
  )
  write_csv(
    unique_fixed_edges,
    file.path(out_dir, "00_unique_fixed_edges_before_fitting.csv")
  )

  all_edges <- all_paths <- all_totals <- list()
  all_boot_edges <- all_boot_paths <- all_boot_totals <- all_boot_fail <- list()
  all_lobo_est <- all_lobo_stab <- list()

  groups <- names(model_spec)
  for (g_idx in seq_along(groups)) {
    group <- groups[g_idx]
    sp <- model_spec[[group]]
    message("\n--- plssem sensitivity: ", analysis_label, ": ", group, " ---")
    d <- prepare_group_data(dat, sp$equations)
    fit <- fit_plssem_group(d, sp$equations)
    edges <- extract_plssem_edges(fit, sp$equations)
    edges$Analysis <- analysis_label
    edges$Group <- group
    edges$Target_label <- unname(target_labels[sp$target])
    edges$Display_significant <- is.finite(edges$P_value) & edges$P_value <= alpha_display
    edges$Fixed_model_path <- TRUE
    edges$Model_specification_basis <- model_spec_basis

    effect_sources <- setdiff(unique(c(edges$From, edges$To)), sp$target)
    paths <- enumerate_path_effects(edges, effect_sources, sp$target)
    paths$Analysis <- analysis_label
    paths$Group <- group
    paths$Target_label <- unname(target_labels[sp$target])
    totals <- total_effects_from_paths(
      paths, sources = effect_sources, target = sp$target
    )
    totals$Analysis <- analysis_label
    totals$Group <- group
    totals$Target_label <- unname(target_labels[sp$target])
    totals$Source_type <- ifelse(
      totals$Source_variable %in% treatments,
      "Treatment", "Intermediate/environmental variable"
    )

    gdir <- file.path(out_dir, "groups", group)
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
    writeLines(fit$syntax, file.path(gdir, "01_plssem_model_syntax.txt"))
    write_csv(edges, file.path(gdir, "02_fixed_path_coefficients.csv"))
    write_csv(paths, file.path(gdir, "03_all_standardized_path_effects.csv"))
    write_csv(totals, file.path(gdir, "04_all_variable_total_effects.csv"))
    saveRDS(fit$fit, file.path(gdir, "plssem_fit.rds"))
    writeLines(capture.output(print(summary(fit$fit))), file.path(gdir, "plssem_summary.txt"))
    write_csv(
      as.data.frame(plssem::parameter_estimates(fit$fit)),
      file.path(gdir, "05_all_parameter_estimates_including_random_structure.csv")
    )
    safe_write_object(
      tryCatch(plssem::fit_measures(fit$fit), error = function(e) NULL),
      file.path(gdir, "06_fit_measures.csv")
    )

    if (run_bootstrap) {
      message(group, ": design-aware whole-plot bootstrap (", n_boot, " replicates)")
      boot <- bootstrap_plssem_group(
        d, sp$equations, effect_sources, sp$target,
        n_boot = n_boot, seed = seed + g_idx, n_cores = n_cores
      )
      if (!is.null(boot)) {
        if (nrow(boot$edge_summary)) {
          names(boot$edge_summary)[names(boot$edge_summary) == "Boot_mean"] <- "Boot_mean_Std_Estimate"
          boot$edge_summary$Analysis <- analysis_label
          boot$edge_summary$Group <- group
          edges <- merge(
            edges, boot$edge_summary,
            by = c("From", "To", "Analysis", "Group"),
            all.x = TRUE, sort = FALSE
          )
        }
        if (nrow(boot$path_summary)) {
          names(boot$path_summary)[names(boot$path_summary) == "Boot_mean"] <- "Boot_mean_Path_effect"
          boot$path_summary$Analysis <- analysis_label
          boot$path_summary$Group <- group
          paths <- merge(
            paths, boot$path_summary,
            by = c("Source_variable", "Target", "Path", "Analysis", "Group"),
            all.x = TRUE, sort = FALSE
          )
        }
        if (nrow(boot$total_summary)) {
          ci <- boot$total_summary
          names(ci)[names(ci) == "Boot_mean"] <- "Boot_mean_Total_effect"
          names(ci)[names(ci) == "Boot_SE"] <- "Boot_SE_Total_effect"
          names(ci)[names(ci) == "CI_low"] <- "Total_CI_low"
          names(ci)[names(ci) == "CI_high"] <- "Total_CI_high"
          names(ci)[names(ci) == "P_boot"] <- "P_boot_Total_effect"
          ci$Analysis <- analysis_label
          ci$Group <- group
          totals <- merge(
            totals, ci,
            by = c("Source_variable", "Target", "Analysis", "Group"),
            all.x = TRUE, sort = FALSE
          )
        }
        if (nrow(boot$direct_summary)) {
          direct_ci <- boot$direct_summary
          stat_cols <- c("N_boot_valid", "Boot_mean", "Boot_SE", "CI_low", "CI_high", "P_boot")
          names(direct_ci)[match(stat_cols, names(direct_ci))] <- paste0("Direct_", stat_cols)
          direct_ci$Analysis <- analysis_label
          direct_ci$Group <- group
          totals <- merge(
            totals, direct_ci,
            by = c("Source_variable", "Target", "Analysis", "Group"),
            all.x = TRUE, sort = FALSE
          )
        }
        if (nrow(boot$indirect_summary)) {
          indirect_ci <- boot$indirect_summary
          stat_cols <- c("N_boot_valid", "Boot_mean", "Boot_SE", "CI_low", "CI_high", "P_boot")
          names(indirect_ci)[match(stat_cols, names(indirect_ci))] <- paste0("Indirect_", stat_cols)
          indirect_ci$Analysis <- analysis_label
          indirect_ci$Group <- group
          totals <- merge(
            totals, indirect_ci,
            by = c("Source_variable", "Target", "Analysis", "Group"),
            all.x = TRUE, sort = FALSE
          )
        }
        all_boot_edges[[group]] <- boot$edges
        all_boot_paths[[group]] <- boot$paths
        all_boot_totals[[group]] <- boot$totals
        if (nrow(boot$failures)) {
          boot$failures$Group <- group
          all_boot_fail[[group]] <- boot$failures
        }
      }
    }

    if (run_lobo) {
      lobo <- lobo_plssem_group(d, sp$equations, edges)
      if (nrow(lobo$estimates)) {
        lobo$estimates$Analysis <- analysis_label
        lobo$estimates$Group <- group
        all_lobo_est[[group]] <- lobo$estimates
      }
      if (nrow(lobo$stability)) {
        lobo$stability$Analysis <- analysis_label
        lobo$stability$Group <- group
        all_lobo_stab[[group]] <- lobo$stability
      }
    }

    all_edges[[group]] <- edges
    all_paths[[group]] <- paths
    all_totals[[group]] <- totals
  }

  combine_nonempty <- function(x) {
    x <- x[vapply(x, function(z) !is.null(z) && nrow(z) > 0, logical(1))]
    if (!length(x)) data.frame() else do.call(rbind, x)
  }
  edges <- combine_nonempty(all_edges)
  paths <- combine_nonempty(all_paths)
  totals <- combine_nonempty(all_totals)
  combined_dir <- file.path(out_dir, "combined")

  write_csv(edges, file.path(combined_dir, "01_all_fixed_path_coefficients.csv"))
  write_csv(edges, file.path(combined_dir, paste0("02_", figure_label, "b_edges_all_fixed.csv")))
  write_csv(
    edges[edges$Display_significant, , drop = FALSE],
    file.path(combined_dir, paste0("03_", figure_label, "b_edges_significant_only.csv"))
  )
  write_csv(paths, file.path(combined_dir, "04_all_unique_standardized_path_effects.csv"))
  write_csv(totals, file.path(combined_dir, paste0("05_", figure_label, "c_all_variable_total_effects.csv")))
  write_csv(
    totals[totals$Source_variable %in% treatments, , drop = FALSE],
    file.path(combined_dir, paste0("05a_", figure_label, "c_treatment_total_effects.csv"))
  )
  write_csv(
    totals[totals$Source_variable == "warm", , drop = FALSE],
    file.path(combined_dir, paste0("06_", figure_label, "c_warming_total_effects_for_plot.csv"))
  )

  if (run_bootstrap) {
    write_csv(combine_nonempty(all_boot_edges), file.path(out_dir, "bootstrap", "B1_edge_bootstrap_replicates.csv"))
    write_csv(combine_nonempty(all_boot_paths), file.path(out_dir, "bootstrap", "B2_path_bootstrap_replicates.csv"))
    write_csv(combine_nonempty(all_boot_totals), file.path(out_dir, "bootstrap", "B3_total_effect_bootstrap_replicates.csv"))
    write_csv(combine_nonempty(all_boot_fail), file.path(out_dir, "bootstrap", "B4_bootstrap_failures.csv"))
  }
  if (run_lobo) {
    write_csv(combine_nonempty(all_lobo_est), file.path(out_dir, "LOBO_robustness", "S1_LOBO_edge_estimates.csv"))
    write_csv(combine_nonempty(all_lobo_stab), file.path(out_dir, "LOBO_robustness", "S2_LOBO_edge_sign_stability.csv"))
  }

  writeLines(c(
    paste0("Analysis=", analysis_label),
    "Role=sensitivity analysis; not the recommended primary manuscript method",
    paste0("plssem_version=", as.character(utils::packageVersion("plssem"))),
    "Estimator=PLS with observed variables; consistent=FALSE",
    "Path_selection=none within plssem; MODEL_SPEC is fixed before fitting",
    paste0("Model_specification_basis=", model_spec_basis),
    "Primary_manuscript_structure=(1|block/WholePlotID)",
    "plssem_cluster_limit=one clustering variable",
    "Block_adjustment_fitted=K-1 fixed treatment contrasts",
    "Random_structure_fitted=(1|BlockWholePlotID)",
    "BlockWholePlotID=interaction(block,WholePlotID)",
    "Interpretation=closest supported sensitivity model; not an exact random-effects equivalent",
    "Total_effect_sources=all non-target variables in each group MODEL_SPEC",
    "Variables_without_a_directed_path=retained with structural total effect 0",
    "Package_internal_bootstrap=FALSE",
    paste0("Design_aware_cluster_bootstrap=", run_bootstrap),
    paste0("N_boot=", if (run_bootstrap) n_boot else 0L),
    "Bootstrap_unit=whole plots resampled within blocks",
    paste0("LOBO_sensitivity=", run_lobo)
  ), file.path(out_dir, "00_analysis_settings.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "00_sessionInfo.txt"))

  invisible(list(edges = edges, paths = paths, totals = totals))
}
