# ==============================================================================
# Shared engine: design-aware mixed-effects piecewise SEM for the TPL study
#
# Experimental structure represented in every component equation:
#   response ~ fixed parents supplied in MODEL_SPEC + (1 | block/WholePlotID)
#
# This expands to:
#   (1 | block) + (1 | block:WholePlotID)
# and matches the hierarchical random structure used in the main analyses.
#
# This engine deliberately contains NO forward/stepwise selection. The supplied
# MODEL_SPEC is the fitted model, whether it was defined a priori or inherited
# from a previous analysis. Paths hidden from a figure remain in the model.
# ==============================================================================

assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Install required package(s) before running: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

safe_sd <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(out) || out <= 0) NA_real_ else out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

empty_df <- function(...) {
  x <- list(...)
  as.data.frame(x, stringsAsFactors = FALSE)
}

write_csv <- function(x, file) {
  utils::write.csv(x, file, row.names = FALSE, quote = FALSE, na = "")
}

read_and_match_splitplot <- function(data_file, meta_file, treatments) {
  dat <- utils::read.csv(
    data_file, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, row.names = 1
  )
  info <- utils::read.csv(
    meta_file, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_meta <- c(
    "plot", "block", "WholePlotID", "warm", "wetting", "drying", "clip"
  )
  missing_meta <- setdiff(required_meta, names(info))
  if (length(missing_meta)) {
    stop("Metadata missing: ", paste(missing_meta, collapse = ", "))
  }

  idx <- match(rownames(dat), info$plot)
  if (anyNA(idx)) {
    stop(
      "Unmatched plot IDs: ",
      paste(rownames(dat)[is.na(idx)], collapse = ", ")
    )
  }
  info2 <- info[idx, , drop = FALSE]
  stopifnot(identical(as.character(rownames(dat)), as.character(info2$plot)))

  dat$plot <- as.character(info2$plot)
  dat$block <- factor(info2$block)
  dat$WholePlotID <- factor(info2$WholePlotID)
  for (v in treatments) dat[[v]] <- as.numeric(info2[[v]])

  if (any(dat$drying == 1 & dat$wetting == 1, na.rm = TRUE)) {
    stop("Some rows are coded as both drying = 1 and wetting = 1.")
  }

  wholeplot_key <- interaction(
    dat$block, dat$WholePlotID,
    drop = TRUE, lex.order = TRUE
  )
  pair_qc <- do.call(rbind, lapply(
    split(seq_len(nrow(dat)), wholeplot_key),
    function(ii) {
      data.frame(
        WholePlotID = as.character(dat$WholePlotID[ii[1]]),
        block = as.character(dat$block[ii[1]]),
        n = length(ii),
        n_clip_levels = length(unique(dat$clip[ii])),
        stringsAsFactors = FALSE
      )
    }
  ))
  if (any(pair_qc$n != 2L | pair_qc$n_clip_levels != 2L)) {
    stop("Each WholePlotID must contain exactly two subplots with clip 0 and 1.")
  }

  list(data = dat, pair_qc = pair_qc)
}

spec_edges <- function(equations) {
  out <- lapply(names(equations), function(response) {
    parents <- unique(equations[[response]])
    if (!length(parents)) return(NULL)
    data.frame(
      From = parents, To = response,
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty_df(From = character(), To = character()))
  unique(do.call(rbind, out))
}


normalize_correlated_errors <- function(correlated_errors) {
  if (is.null(correlated_errors) || !length(correlated_errors)) return(list())
  if (!is.list(correlated_errors)) {
    stop("correlated_errors must be a list of two-variable character vectors.")
  }

  out <- lapply(seq_along(correlated_errors), function(i) {
    pair <- as.character(correlated_errors[[i]])
    if (length(pair) != 2L || anyNA(pair) || any(!nzchar(pair))) {
      stop(
        "Each correlated error must contain exactly two non-empty variable names; ",
        "problem at correlated_errors[[", i, "]]."
      )
    }
    if (identical(pair[1], pair[2])) {
      stop("A correlated error cannot pair a variable with itself: ", pair[1])
    }
    pair
  })

  # Treat X ~~ Y and Y ~~ X as the same specification.
  keys <- vapply(
    out,
    function(pair) paste(sort(pair), collapse = "||"),
    character(1)
  )
  out[!duplicated(keys)]
}

correlated_error_vars <- function(correlated_errors) {
  ce <- normalize_correlated_errors(correlated_errors)
  if (!length(ce)) return(character())
  unique(unlist(ce, use.names = FALSE))
}

make_correlated_error_terms <- function(correlated_errors) {
  ce <- normalize_correlated_errors(correlated_errors)
  lapply(ce, function(pair) {
    # piecewiseSEM::`%~~%` returns a character string with class
    # "formula.cerror". Constructing the same object directly is robust inside
    # do.call() and avoids non-standard-evaluation problems.
    term <- paste(pair[1], "~~", pair[2])
    class(term) <- "formula.cerror"
    term
  })
}

spec_correlated_errors <- function(model_spec) {
  rows <- list()
  k <- 1L
  for (group in names(model_spec)) {
    ce <- normalize_correlated_errors(model_spec[[group]]$correlated_errors %||% NULL)
    if (!length(ce)) next
    for (pair in ce) {
      rows[[k]] <- data.frame(
        Group = group,
        Variable1 = pair[1],
        Variable2 = pair[2],
        Correlated_error = paste(pair[1], "~~", pair[2]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (!length(rows)) {
    return(empty_df(
      Group = character(),
      Variable1 = character(),
      Variable2 = character(),
      Correlated_error = character()
    ))
  }
  do.call(rbind, rows)
}

validate_model_spec <- function(model_spec, data_names, treatments) {
  if (!length(model_spec)) stop("MODEL_SPEC is empty.")

  for (group in names(model_spec)) {
    sp <- model_spec[[group]]
    if (is.null(sp$target) || is.null(sp$equations)) {
      stop(group, ": each specification needs target and equations.")
    }
    if (!sp$target %in% names(sp$equations)) {
      stop(group, ": target is not an endogenous response in equations.")
    }

    edges <- spec_edges(sp$equations)
    cerrors <- normalize_correlated_errors(sp$correlated_errors %||% NULL)
    cerror_variables <- correlated_error_vars(cerrors)

    needed <- unique(c(
      edges$From, edges$To, cerror_variables,
      "block", "WholePlotID", "clip", "plot"
    ))
    missing <- setdiff(needed, data_names)
    if (length(missing)) {
      stop(group, ": data missing variable(s): ", paste(missing, collapse = ", "))
    }
    if (any(edges$From == edges$To)) stop(group, ": self-loop in MODEL_SPEC.")

    # Kahn topological-sort check applies only to directed paths. Correlated
    # errors are non-directional and therefore do not create DAG cycles.
    nodes <- unique(c(edges$From, edges$To))
    indegree <- setNames(integer(length(nodes)), nodes)
    for (v in edges$To) indegree[v] <- indegree[v] + 1L
    remaining <- edges
    queue <- names(indegree)[indegree == 0L]
    visited <- character()
    while (length(queue)) {
      node <- queue[1]
      queue <- queue[-1]
      visited <- c(visited, node)
      children <- remaining$To[remaining$From == node]
      remaining <- remaining[remaining$From != node, , drop = FALSE]
      for (child in children) {
        indegree[child] <- indegree[child] - 1L
        if (indegree[child] == 0L) queue <- c(queue, child)
      }
    }
    if (length(unique(visited)) != length(nodes)) {
      stop(group, ": MODEL_SPEC directed paths must form a directed acyclic graph.")
    }

    if (!any(treatments %in% edges$From)) {
      warning(group, ": no treatment appears as an exogenous path source.")
    }
  }
  invisible(TRUE)
}

prepare_group_data <- function(
    dat,
    equations,
    correlated_errors = NULL,
    min_wholeplots_per_block = 3L) {

  edges <- spec_edges(equations)
  cerror_variables <- correlated_error_vars(correlated_errors)
  vars <- unique(c(
    edges$From, edges$To, cerror_variables,
    "plot", "block", "WholePlotID", "clip"
  ))
  d <- dat[, vars, drop = FALSE]

  numeric_vars <- setdiff(vars, c("plot", "block", "WholePlotID"))
  for (v in numeric_vars) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[stats::complete.cases(d), , drop = FALSE]

  # If one member of a whole-plot pair is incomplete, discard the entire pair.
  pair_key <- as.character(interaction(
    d$block, d$WholePlotID,
    drop = TRUE, lex.order = TRUE
  ))
  good <- vapply(split(seq_len(nrow(d)), pair_key), function(ii) {
    length(ii) == 2L && length(unique(d$clip[ii])) == 2L
  }, logical(1))
  d <- d[pair_key %in% names(good)[good], , drop = FALSE]
  d$block <- droplevels(factor(d$block))
  d$WholePlotID <- droplevels(factor(d$WholePlotID))

  # Recompute after filtering so the key is aligned with the retained rows.
  pair_key <- as.character(interaction(
    d$block, d$WholePlotID,
    drop = TRUE, lex.order = TRUE
  ))
  wp_by_block <- tapply(pair_key, d$block, function(x) length(unique(x)))
  if (length(wp_by_block) < 3L || any(wp_by_block < min_wholeplots_per_block)) {
    stop(
      "Too few complete whole plots after complete-case filtering: ",
      paste(names(wp_by_block), wp_by_block, sep = "=", collapse = ", ")
    )
  }
  d
}

make_component_formula <- function(response, parents, env = parent.frame()) {
  rhs <- paste(
    c(unique(parents), "(1 | block/WholePlotID)"),
    collapse = " + "
  )
  stats::as.formula(paste(response, "~", rhs), env = env)
}

fit_component_models <- function(dat, equations) {
  out <- lapply(names(equations), function(response) {
    # piecewiseSEM's d-separation tests update() each lmer model with
    # `data = data`. update.merMod evaluates that call in the formula
    # environment. If the formula was created inside this function but that
    # environment contains only an object called `dat`, R can instead resolve
    # `data` to utils::data(), producing:
    #   list2env(data): first argument must be a named list
    # Give every component model its own persistent formula environment with
    # both names bound to the complete analysis data. This also keeps saved RDS
    # fits reproducible after fit_component_models() returns.
    model_env <- new.env(parent = parent.frame())
    model_env$data <- dat
    model_env$dat <- dat
    f <- make_component_formula(
      response, equations[[response]], env = model_env
    )
    model_env$f <- f

    # Evaluate inside model_env so the stored call uses stable `f` and `data`
    # bindings. REML is required for the F/Kenward-Roger tests used by
    # piecewiseSEM::coefs() for lme4 models.
    eval(
      quote(lme4::lmer(
        f,
        data = data,
        REML = TRUE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5),
          check.conv.singular = "ignore"
        )
      )),
      envir = model_env
    )
  })
  names(out) <- names(equations)
  out
}

extract_edge_estimates <- function(models, equations, dat) {
  rows <- list()
  k <- 1L
  for (response in names(equations)) {
    sm <- coef(summary(models[[response]]))
    sy <- safe_sd(dat[[response]])
    for (parent in equations[[response]]) {
      if (!parent %in% rownames(sm)) {
        stop("Coefficient absent from component model: ", parent, " -> ", response)
      }
      sx <- safe_sd(dat[[parent]])
      scale_ratio <- sx / sy
      estimate <- unname(sm[parent, "Estimate"])
      std_error <- unname(sm[parent, "Std. Error"])
      z_value <- estimate / std_error
      rows[[k]] <- data.frame(
        From = parent,
        To = response,
        Estimate_raw = estimate,
        SE_raw = std_error,
        Std_Estimate = estimate * scale_ratio,
        Std_SE = std_error * scale_ratio,
        P_value_normal_approx = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

extract_piecewise_pvalues <- function(sem) {
  tab <- tryCatch(
    as.data.frame(piecewiseSEM::coefs(
      sem, standardize = "scale", test.statistic = "F", test.type = "II"
    )),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab)) return(NULL)
  needed <- c("Response", "Predictor", "P.Value")
  if (!all(needed %in% names(tab))) return(NULL)
  out <- tab[, needed, drop = FALSE]
  names(out) <- c("To", "From", "P_value_KR")
  out
}


safe_numeric_column <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

extract_correlated_error_estimates <- function(sem, correlated_errors) {
  ce <- normalize_correlated_errors(correlated_errors)
  if (!length(ce)) {
    return(empty_df(
      From = character(), To = character(), Estimate_raw = numeric(),
      Std_Estimate = numeric(), P_value = numeric(), P_method = character(),
      Direction = character()
    ))
  }

  tab <- tryCatch(
    as.data.frame(piecewiseSEM::coefs(
      sem, standardize = "scale", test.statistic = "F", test.type = "II"
    )),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab)) {
    stop("Could not extract piecewiseSEM coefficients for correlated errors.")
  }
  needed <- c("Response", "Predictor", "Estimate", "P.Value", "Std.Estimate")
  if (!all(needed %in% names(tab))) {
    stop(
      "piecewiseSEM::coefs() is missing correlated-error columns: ",
      paste(setdiff(needed, names(tab)), collapse = ", ")
    )
  }

  response_clean <- sub("^~~", "", as.character(tab$Response))
  predictor_clean <- sub("^~~", "", as.character(tab$Predictor))

  rows <- lapply(ce, function(pair) {
    idx <- which(
      (response_clean == pair[1] & predictor_clean == pair[2]) |
        (response_clean == pair[2] & predictor_clean == pair[1])
    )
    if (length(idx) != 1L) {
      stop(
        "Expected exactly one coefficient row for correlated error ",
        paste(pair, collapse = " ~~ "), "; found ", length(idx), "."
      )
    }
    i <- idx[1]
    est <- safe_numeric_column(tab$Estimate[i])
    std_est <- safe_numeric_column(tab$Std.Estimate[i])
    pval <- safe_numeric_column(tab$P.Value[i])
    data.frame(
      From = pair[1],
      To = pair[2],
      Estimate_raw = est,
      Std_Estimate = std_est,
      P_value = pval,
      P_method = "piecewiseSEM correlated error",
      Direction = if (!is.finite(std_est)) NA_character_ else if (std_est >= 0) "positive" else "negative",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

significance_code <- function(p) {
  out <- rep("", length(p))
  good <- is.finite(p)
  out[good & p < 0.001] <- "***"
  out[good & p >= 0.001 & p < 0.01] <- "**"
  out[good & p >= 0.01 & p < 0.05] <- "*"
  out
}

pretty_TPL_variable_label <- function(x) {
  labels <- c(
    warm = "Warming",
    drying = "Reduced precipitation",
    wetting = "Increased precipitation",
    clip = "Clipping",
    temperature_mean = "Temperature mean",
    moisture_mean = "Soil moisture mean",
    moisture_cv = "Soil moisture CV",
    pH_mean = "Soil pH mean",
    NO3N_mean = "NO3-N mean",
    NO3N_cv = "NO3-N CV",
    TC_mean = "Total C mean",
    TC_cv = "Total C CV",
    plant_biomass_mean = "Plant biomass mean",
    plant_biomass_cv = "Plant biomass CV",
    plant_richness_mean = "Plant richness mean",
    plant_richness_cv = "Plant richness CV"
  )
  out <- as.character(x)
  hit <- out %in% names(labels)
  out[hit] <- unname(labels[out[hit]])
  out
}

column_or <- function(x, name, default = NA_real_) {
  if (name %in% names(x)) return(x[[name]])
  rep(default, nrow(x))
}

make_path_diagram_plot_table <- function(
    edges,
    correlated_error_edges = NULL,
    lobo_stability = NULL,
    alpha_display = 0.05,
    figure_label = "Figure") {

  if (is.null(edges) || !nrow(edges)) return(data.frame())

  directed <- data.frame(
    Analysis = edges$Analysis,
    Figure = figure_label,
    Group = edges$Group,
    Target_label = edges$Target_label,
    Edge_type = "directed",
    From = edges$From,
    To = edges$To,
    Std_Estimate = edges$Std_Estimate,
    P_value = edges$P_value,
    P_method = edges$P_method,
    Direction = edges$Direction,
    Significant = is.finite(edges$P_value) & edges$P_value <= alpha_display,
    Significance = significance_code(edges$P_value),
    N_boot_valid = column_or(edges, "N_boot_valid", NA_integer_),
    Boot_mean_Std_Estimate = column_or(edges, "Boot_mean_Std_Estimate"),
    Boot_CI_low = column_or(edges, "CI_low"),
    Boot_CI_high = column_or(edges, "CI_high"),
    P_boot = column_or(edges, "P_boot"),
    stringsAsFactors = FALSE
  )
  directed$Bootstrap_supported <- with(
    directed,
    is.finite(Boot_CI_low) & is.finite(Boot_CI_high) &
      ((Boot_CI_low > 0 & Boot_CI_high > 0) |
         (Boot_CI_low < 0 & Boot_CI_high < 0))
  )

  if (!is.null(lobo_stability) && nrow(lobo_stability)) {
    keep <- intersect(
      c(
        "From", "To", "Analysis", "Group", "LOBO_min", "LOBO_max",
        "LOBO_mean", "Sign_consistency_vs_full", "N_LOBO_folds"
      ),
      names(lobo_stability)
    )
    directed <- merge(
      directed,
      lobo_stability[, keep, drop = FALSE],
      by = c("From", "To", "Analysis", "Group"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    directed$LOBO_min <- NA_real_
    directed$LOBO_max <- NA_real_
    directed$LOBO_mean <- NA_real_
    directed$Sign_consistency_vs_full <- NA_real_
    directed$N_LOBO_folds <- NA_integer_
  }

  cerror_plot <- data.frame()
  if (!is.null(correlated_error_edges) && nrow(correlated_error_edges)) {
    cerror_plot <- data.frame(
      Analysis = correlated_error_edges$Analysis,
      Figure = figure_label,
      Group = correlated_error_edges$Group,
      Target_label = correlated_error_edges$Target_label,
      Edge_type = "correlated_error",
      From = correlated_error_edges$From,
      To = correlated_error_edges$To,
      Std_Estimate = correlated_error_edges$Std_Estimate,
      P_value = correlated_error_edges$P_value,
      P_method = correlated_error_edges$P_method,
      Direction = correlated_error_edges$Direction,
      Significant = is.finite(correlated_error_edges$P_value) &
        correlated_error_edges$P_value <= alpha_display,
      Significance = significance_code(correlated_error_edges$P_value),
      N_boot_valid = NA_integer_,
      Boot_mean_Std_Estimate = NA_real_,
      Boot_CI_low = NA_real_,
      Boot_CI_high = NA_real_,
      P_boot = NA_real_,
      Bootstrap_supported = NA,
      LOBO_min = NA_real_,
      LOBO_max = NA_real_,
      LOBO_mean = NA_real_,
      Sign_consistency_vs_full = NA_real_,
      N_LOBO_folds = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  out <- if (nrow(cerror_plot)) rbind(directed, cerror_plot) else directed
  out$From_label <- pretty_TPL_variable_label(out$From)
  out$To_label <- pretty_TPL_variable_label(out$To)
  target_rows <- grepl("_beta$", out$To)
  out$To_label[target_rows] <- out$Target_label[target_rows]
  out$Arrow <- out$Edge_type == "directed"
  out$Line_type <- ifelse(out$Edge_type == "correlated_error", "dashed", "solid")
  out$Plot_include <- out$Significant
  out$Edge_label <- ifelse(
    is.finite(out$Std_Estimate),
    paste0(sprintf("%.2f", out$Std_Estimate), out$Significance),
    ""
  )

  preferred <- c(
    "Analysis", "Figure", "Group", "Target_label", "Edge_type",
    "From", "From_label", "To", "To_label", "Std_Estimate",
    "P_value", "P_method", "Significant", "Significance", "Direction",
    "N_boot_valid", "Boot_mean_Std_Estimate", "Boot_CI_low", "Boot_CI_high",
    "P_boot", "Bootstrap_supported", "LOBO_min", "LOBO_max", "LOBO_mean",
    "Sign_consistency_vs_full", "N_LOBO_folds", "Arrow", "Line_type",
    "Plot_include", "Edge_label"
  )
  out <- out[, preferred, drop = FALSE]
  out[order(out$Group, out$Edge_type, out$From, out$To), , drop = FALSE]
}

make_effect_plot_table <- function(
    totals,
    plot_source_order = NULL,
    figure_label = "Figure") {

  if (is.null(totals) || !nrow(totals)) return(data.frame())

  if (is.null(plot_source_order)) {
    plot_source_order <- unique(as.character(totals$Source_variable))
  }
  extra <- setdiff(unique(as.character(totals$Source_variable)), plot_source_order)
  plot_source_order <- unique(c(plot_source_order, extra))

  out <- data.frame(
    Analysis = totals$Analysis,
    Figure = figure_label,
    Group = totals$Group,
    Target_label = totals$Target_label,
    Source_variable = totals$Source_variable,
    Source_label = pretty_TPL_variable_label(totals$Source_variable),
    Source_type = totals$Source_type,
    Source_order = match(totals$Source_variable, plot_source_order),
    Direct_effect = totals$Direct_effect,
    Indirect_effect = totals$Indirect_effect,
    Total_effect = totals$Total_effect,
    N_unique_paths = totals$N_unique_paths,
    Has_directed_path = totals$Has_directed_path,
    stringsAsFactors = FALSE
  )

  out$Direct_boot_mean <- column_or(totals, "Direct_Boot_mean")
  out$Direct_CI_low <- column_or(totals, "Direct_CI_low")
  out$Direct_CI_high <- column_or(totals, "Direct_CI_high")
  out$Direct_P_boot <- column_or(totals, "Direct_P_boot")
  out$Indirect_boot_mean <- column_or(totals, "Indirect_Boot_mean")
  out$Indirect_CI_low <- column_or(totals, "Indirect_CI_low")
  out$Indirect_CI_high <- column_or(totals, "Indirect_CI_high")
  out$Indirect_P_boot <- column_or(totals, "Indirect_P_boot")
  out$Total_boot_mean <- column_or(totals, "Boot_mean_Total_effect")
  out$Total_boot_SE <- column_or(totals, "Boot_SE_Total_effect")
  out$Total_CI_low <- column_or(totals, "Total_CI_low")
  out$Total_CI_high <- column_or(totals, "Total_CI_high")
  out$Total_P_boot <- column_or(totals, "P_boot_Total_effect")
  out$N_boot_valid <- column_or(totals, "N_boot_valid", NA_integer_)

  out$Total_bootstrap_supported <- with(
    out,
    is.finite(Total_CI_low) & is.finite(Total_CI_high) &
      ((Total_CI_low > 0 & Total_CI_high > 0) |
         (Total_CI_low < 0 & Total_CI_high < 0))
  )
  out$Total_direction <- ifelse(
    out$Total_effect > 0, "positive",
    ifelse(out$Total_effect < 0, "negative", "zero")
  )
  # User-requested plotting table keeps every treatment and every intermediate
  # environmental variable, including structural zero effects, so that panels
  # can compare the same source set across taxa.
  out$Plot_include <- TRUE

  out[order(out$Group, out$Source_order, out$Source_variable), , drop = FALSE]
}

component_diagnostics <- function(models, equations, dat) {
  do.call(rbind, lapply(names(models), function(response) {
    m <- models[[response]]
    conv <- m@optinfo$conv$lme4$messages
    vc <- as.data.frame(lme4::VarCorr(m))
    grp <- as.character(vc$grp)
    intercept_var <- is.na(vc$var2)
    block_var <- vc$vcov[grp == "block" & intercept_var]
    block_wp_var <- vc$vcov[
      grepl("block", grp, fixed = TRUE) &
        grepl("WholePlotID", grp, fixed = TRUE) &
        intercept_var
    ]
    r2 <- tryCatch(performance::r2_nakagawa(m), error = function(e) NULL)
    data.frame(
      Response = response,
      N_subplots = stats::nobs(m),
      N_wholeplots = length(unique(interaction(
        dat$block, dat$WholePlotID,
        drop = TRUE, lex.order = TRUE
      ))),
      N_blocks = length(unique(dat$block)),
      N_substantive_predictors = length(equations[[response]]),
      Substantive_predictors = paste(equations[[response]], collapse = ";"),
      Singular = lme4::isSingular(m, tol = 1e-5),
      Block_variance = if (length(block_var)) block_var[1] else NA_real_,
      Block_WholePlot_variance = if (length(block_wp_var)) block_wp_var[1] else NA_real_,
      R2_marginal = if (is.null(r2)) NA_real_ else unname(r2$R2_marginal),
      R2_conditional = if (is.null(r2)) NA_real_ else unname(r2$R2_conditional),
      Convergence_message = if (is.null(conv)) "" else paste(conv, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
}

fit_piecewise_group <- function(dat, equations, correlated_errors = NULL) {
  models <- fit_component_models(dat, equations)
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("Component-model list must retain response-variable names.")
  }

  cerror_terms <- make_correlated_error_terms(correlated_errors)

  sem <- do.call(
    piecewiseSEM::psem,
    # Retain response names for downstream bookkeeping. psem() also receives
    # the complete data frame explicitly; component calls retain their own
    # stable formula/data environments for update-based d-separation tests.
    c(models, cerror_terms, list(data = dat))
  )

  # Fail loudly if a requested correlated error was silently dropped. This
  # protects the d-separation basis set and Fisher's C from being computed on
  # a model different from MODEL_SPEC.
  if (length(cerror_terms)) {
    requested <- gsub(
      "[[:space:]]+", "",
      vapply(cerror_terms, as.character, character(1))
    )
    present_idx <- vapply(
      sem,
      function(x) inherits(x, "formula.cerror"),
      logical(1)
    )
    present <- gsub(
      "[[:space:]]+", "",
      vapply(sem[present_idx], as.character, character(1))
    )
    missing_cerrors <- setdiff(requested, present)
    if (length(missing_cerrors)) {
      stop(
        "Requested correlated error(s) were not retained by piecewiseSEM::psem(): ",
        paste(missing_cerrors, collapse = ", ")
      )
    }

    # Stronger functional check: a correlated-error pair must be removed from
    # the directed-separation basis set. This catches cases where the object is
    # present in the psem list but is not actually being honored by basisSet().
    basis <- piecewiseSEM::basisSet(sem)
    ce_pairs <- normalize_correlated_errors(correlated_errors)
    still_tested <- vapply(ce_pairs, function(pair) {
      any(vapply(basis, function(claim) {
        length(claim) >= 2L && all(pair %in% claim[1:2])
      }, logical(1)))
    }, logical(1))
    if (any(still_tested)) {
      bad <- vapply(
        ce_pairs[still_tested],
        function(pair) paste(pair, collapse = " ~~ "),
        character(1)
      )
      stop(
        "Correlated error(s) remain in the d-separation basis set: ",
        paste(bad, collapse = ", "),
        ". Check the installed piecewiseSEM version."
      )
    }
  }

  edges <- extract_edge_estimates(models, equations, dat)
  ptab <- extract_piecewise_pvalues(sem)
  if (!is.null(ptab)) {
    edges <- merge(edges, ptab, by = c("From", "To"), all.x = TRUE, sort = FALSE)
  } else {
    edges$P_value_KR <- NA_real_
  }
  edges$P_value <- ifelse(
    is.finite(edges$P_value_KR),
    edges$P_value_KR,
    edges$P_value_normal_approx
  )
  edges$P_method <- ifelse(
    is.finite(edges$P_value_KR),
    "piecewiseSEM Kenward-Roger",
    "normal approximation fallback"
  )
  edges$Direction <- ifelse(edges$Std_Estimate >= 0, "positive", "negative")
  cerror_edges <- extract_correlated_error_estimates(sem, correlated_errors)

  list(
    models = models,
    sem = sem,
    edges = edges,
    cerror_edges = cerror_edges,
    correlated_errors = normalize_correlated_errors(correlated_errors),
    diagnostics = component_diagnostics(models, equations, dat)
  )
}

find_simple_paths <- function(edges, source, target, max_depth = 8L) {
  if (!nrow(edges)) return(list())
  adjacency <- split(edges$To, edges$From)
  out <- list()
  dfs <- function(node, path) {
    if (length(path) - 1L > max_depth) return(invisible(NULL))
    if (identical(node, target)) {
      out[[length(out) + 1L]] <<- path
      return(invisible(NULL))
    }
    next_nodes <- adjacency[[node]]
    if (is.null(next_nodes)) return(invisible(NULL))
    for (next_node in unique(next_nodes)) {
      if (next_node %in% path) next
      dfs(next_node, c(path, next_node))
    }
    invisible(NULL)
  }
  dfs(source, source)
  out
}

enumerate_path_effects <- function(edges, sources, target, estimate_col = "Std_Estimate") {
  records <- list()
  k <- 1L
  for (source in sources) {
    paths <- find_simple_paths(edges, source, target)
    if (!length(paths)) next
    for (path in paths) {
      vals <- numeric(length(path) - 1L)
      for (j in seq_len(length(path) - 1L)) {
        ee <- edges[edges$From == path[j] & edges$To == path[j + 1L], , drop = FALSE]
        if (nrow(ee) != 1L) {
          stop("Expected exactly one edge: ", path[j], " -> ", path[j + 1L])
        }
        vals[j] <- ee[[estimate_col]][1]
      }
      n_edges <- length(vals)
      records[[k]] <- data.frame(
        Source_variable = source,
        Target = target,
        Path = paste(path, collapse = " -> "),
        N_edges = n_edges,
        Path_type = ifelse(n_edges == 1L, "Direct", "Indirect"),
        Path_effect = prod(vals),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (!length(records)) {
    return(empty_df(
      Source_variable = character(), Target = character(), Path = character(),
      N_edges = integer(), Path_type = character(), Path_effect = numeric()
    ))
  }
  do.call(rbind, records)
}

total_effects_from_paths <- function(paths, sources = NULL, target = NULL) {
  if (!is.null(sources) && !is.null(target)) {
    keys <- data.frame(
      Source_variable = unique(as.character(sources)),
      Target = as.character(target),
      stringsAsFactors = FALSE
    )
  } else if (nrow(paths)) {
    keys <- unique(paths[, c("Source_variable", "Target"), drop = FALSE])
  } else {
    return(empty_df(
      Source_variable = character(), Target = character(),
      Direct_effect = numeric(), Indirect_effect = numeric(), Total_effect = numeric(),
      N_unique_paths = integer(), Has_directed_path = logical()
    ))
  }
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    dd <- paths[
      paths$Source_variable == keys$Source_variable[i] &
        paths$Target == keys$Target[i], , drop = FALSE
    ]
    direct <- sum(dd$Path_effect[dd$Path_type == "Direct"], na.rm = TRUE)
    indirect <- sum(dd$Path_effect[dd$Path_type == "Indirect"], na.rm = TRUE)
    data.frame(
      Source_variable = keys$Source_variable[i],
      Target = keys$Target[i],
      Direct_effect = direct,
      Indirect_effect = indirect,
      Total_effect = direct + indirect,
      N_unique_paths = nrow(dd),
      Has_directed_path = nrow(dd) > 0L,
      stringsAsFactors = FALSE
    )
  }))
}

resample_wholeplots_within_blocks <- function(dat) {
  pieces <- list()
  k <- 1L
  block_levels <- levels(dat$block)
  for (b in block_levels) {
    db <- dat[dat$block == b, , drop = FALSE]
    wps <- unique(as.character(db$WholePlotID))
    donors <- sample(wps, length(wps), replace = TRUE)
    for (j in seq_along(donors)) {
      x <- db[as.character(db$WholePlotID) == donors[j], , drop = FALSE]
      x$WholePlotID <- paste0("boot_", b, "_", j)
      pieces[[k]] <- x
      k <- k + 1L
    }
  }
  out <- do.call(rbind, pieces)
  out$block <- factor(out$block, levels = block_levels)
  out$WholePlotID <- factor(out$WholePlotID)
  rownames(out) <- NULL
  out
}

summarize_bootstrap <- function(dat, key_cols, value_col) {
  if (is.null(dat) || !nrow(dat)) return(data.frame())
  key <- interaction(dat[, key_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  pieces <- split(dat, key)
  do.call(rbind, lapply(pieces, function(dd) {
    x <- dd[[value_col]]
    x <- x[is.finite(x)]
    base <- dd[1, key_cols, drop = FALSE]
    if (!length(x)) {
      base$N_boot_valid <- 0L
      base$Boot_mean <- base$Boot_SE <- base$CI_low <- base$CI_high <- base$P_boot <- NA_real_
      return(base)
    }
    p_lower <- (1 + sum(x <= 0)) / (length(x) + 1)
    p_upper <- (1 + sum(x >= 0)) / (length(x) + 1)
    base$N_boot_valid <- length(x)
    base$Boot_mean <- mean(x)
    base$Boot_SE <- stats::sd(x)
    ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 6)
    base$CI_low <- ci[1]
    base$CI_high <- ci[2]
    base$P_boot <- min(1, 2 * min(p_lower, p_upper))
    base
  }))
}

bootstrap_group <- function(dat, equations, sources, target, n_boot, seed, n_cores) {
  if (n_boot <= 0L) return(NULL)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_boot)

  one <- function(i) {
    set.seed(seeds[i])
    dboot <- resample_wholeplots_within_blocks(dat)
    ans <- tryCatch({
      models <- fit_component_models(dboot, equations)
      edges <- extract_edge_estimates(models, equations, dboot)
      paths <- enumerate_path_effects(edges, sources, target)
      totals <- total_effects_from_paths(paths, sources = sources, target = target)
      edges$Replicate <- i
      paths$Replicate <- i
      totals$Replicate <- i
      list(edges = edges, paths = paths, totals = totals, error = NULL)
    }, error = function(e) {
      list(edges = NULL, paths = NULL, totals = NULL, error = conditionMessage(e))
    })
    ans
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

lobo_group <- function(dat, equations, full_edges) {
  blocks <- levels(dat$block)
  rows <- list()
  k <- 1L
  for (b in blocks) {
    dtrain <- droplevels(dat[dat$block != b, , drop = FALSE])
    ans <- tryCatch({
      models <- fit_component_models(dtrain, equations)
      extract_edge_estimates(models, equations, dtrain)
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

safe_write_object <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.atomic(x) && is.null(dim(x))) x <- data.frame(Value = x)
  x <- tryCatch(as.data.frame(x), error = function(e) data.frame(Text = capture.output(print(x))))
  write_csv(x, file)
  invisible(TRUE)
}

capture_psem_summary <- function(sem, models) {
  main <- tryCatch(
    capture.output(print(summary(sem, conserve = TRUE))),
    error = function(e) e
  )
  if (!inherits(main, "error")) return(main)

  capture_section <- function(title, x) {
    body <- tryCatch(
      capture.output(print(force(x))),
      error = function(e) paste("FAILED:", conditionMessage(e))
    )
    c(paste0("----- ", title, " -----"), body, "")
  }

  out <- c(
    paste("summary.psem failed:", conditionMessage(main)),
    "A component-level fallback report follows so this file is still diagnostic.",
    ""
  )
  out <- c(
    out,
    capture_section(
      "piecewiseSEM coefficients",
      piecewiseSEM::coefs(
        sem, standardize = "scale", test.statistic = "F", test.type = "II"
      )
    ),
    capture_section("Fisher's C", piecewiseSEM::fisherC(sem, conserve = TRUE)),
    capture_section("d-separation claims", piecewiseSEM::dSep(sem, conserve = TRUE)),
    capture_section("d-separation AIC", stats::AIC(sem, AIC.type = "dsep"))
  )
  for (response in names(models)) {
    out <- c(
      out,
      capture_section(
        paste0("component model: ", response),
        summary(models[[response]])
      )
    )
  }
  out
}

run_splitplot_piecewise_sem <- function(
    data_file,
    meta_file,
    out_dir,
    analysis_label,
    figure_label,
    model_spec,
    target_labels,
    model_spec_basis = "fixed paths supplied before piecewiseSEM fitting",
    treatments = c("warm", "drying", "wetting", "clip"),
    alpha_display = 0.05,
    run_bootstrap = TRUE,
    n_boot = 999L,
    seed = 8302026L,
    n_cores = 1L,
    run_lobo = TRUE,
    plot_source_order = NULL) {

  assert_packages(c("lme4", "piecewiseSEM", "performance"))
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
  correlated_errors_by_group <- spec_correlated_errors(model_spec)

  design_qc <- data.frame(
    Analysis = analysis_label,
    N_subplots = nrow(dat),
    N_blocks = length(unique(dat$block)),
    N_wholeplots = length(unique(interaction(
      dat$block, dat$WholePlotID,
      drop = TRUE, lex.order = TRUE
    ))),
    Block_treatment_unit = "WholePlotID",
    Clipping_treatment_unit = "subplot within WholePlotID",
    Component_model_structure = "fixed parents + (1|block/WholePlotID)",
    Model_specification_basis = model_spec_basis,
    N_unique_fixed_edges = nrow(unique_fixed_edges),
    N_correlated_errors = nrow(correlated_errors_by_group),
    stringsAsFactors = FALSE
  )
  write_csv(design_qc, file.path(out_dir, "00_design_QC.csv"))
  write_csv(input$pair_qc, file.path(out_dir, "00_wholeplot_pair_QC.csv"))
  write_csv(
    fixed_paths,
    file.path(out_dir, "00_fixed_paths_by_group_before_fitting.csv")
  )
  write_csv(
    unique_fixed_edges,
    file.path(out_dir, "00_unique_fixed_edges_before_fitting.csv")
  )
  write_csv(
    correlated_errors_by_group,
    file.path(out_dir, "00_correlated_errors_by_group_before_fitting.csv")
  )

  all_edges <- list()
  all_cerror_edges <- list()
  all_paths <- list()
  all_totals <- list()
  all_diag <- list()
  all_fit <- list()
  all_boot_edges <- list()
  all_boot_paths <- list()
  all_boot_totals <- list()
  all_boot_fail <- list()
  all_lobo_est <- list()
  all_lobo_stab <- list()

  group_names <- names(model_spec)
  for (g_idx in seq_along(group_names)) {
    group <- group_names[g_idx]
    message("\n--- ", analysis_label, ": ", group, " ---")
    sp <- model_spec[[group]]
    cerrors <- sp$correlated_errors %||% NULL
    d <- prepare_group_data(
      dat,
      sp$equations,
      correlated_errors = cerrors
    )
    fit <- fit_piecewise_group(
      d,
      sp$equations,
      correlated_errors = cerrors
    )

    edges <- fit$edges
    edges$Analysis <- analysis_label
    edges$Group <- group
    edges$Target_label <- unname(target_labels[sp$target])
    edges$Display_significant <- is.finite(edges$P_value) & edges$P_value <= alpha_display
    edges$Fixed_model_path <- TRUE
    edges$Model_specification_basis <- model_spec_basis

    cerror_edges <- fit$cerror_edges
    if (!is.null(cerror_edges) && nrow(cerror_edges)) {
      cerror_edges$Analysis <- analysis_label
      cerror_edges$Group <- group
      cerror_edges$Target_label <- unname(target_labels[sp$target])
      cerror_edges$Display_significant <-
        is.finite(cerror_edges$P_value) & cerror_edges$P_value <= alpha_display
      cerror_edges$Model_specification_basis <- model_spec_basis
    }

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

    diag <- fit$diagnostics
    diag$Analysis <- analysis_label
    diag$Group <- group

    gdir <- file.path(out_dir, "groups", group)
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
    write_csv(edges, file.path(gdir, "01_fixed_path_coefficients.csv"))
    write_csv(
      cerror_edges,
      file.path(gdir, "01a_correlated_error_coefficients.csv")
    )
    write_csv(diag, file.path(gdir, "02_component_model_diagnostics_R2.csv"))
    write_csv(paths, file.path(gdir, "03_all_standardized_path_effects.csv"))
    write_csv(totals, file.path(gdir, "04_all_variable_total_effects.csv"))
    group_cerrors <- correlated_errors_by_group[
      correlated_errors_by_group$Group == group, , drop = FALSE
    ]
    write_csv(group_cerrors, file.path(gdir, "00_correlated_errors_specified.csv"))
    saveRDS(fit$sem, file.path(gdir, "piecewiseSEM_fit.rds"))
    summary_text <- capture_psem_summary(fit$sem, fit$models)
    writeLines(summary_text, file.path(gdir, "piecewiseSEM_summary.txt"))
    safe_write_object(
      tryCatch(piecewiseSEM::fisherC(fit$sem, conserve = TRUE), error = function(e) NULL),
      file.path(gdir, "05_Fisher_C_global_fit.csv")
    )
    safe_write_object(
      tryCatch(stats::AIC(fit$sem, AIC.type = "dsep"), error = function(e) NULL),
      file.path(gdir, "06_dsep_AIC.csv")
    )

    if (run_bootstrap) {
      message(group, ": whole-plot cluster bootstrap (", n_boot, " replicates)")
      boot <- bootstrap_group(
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
          total_ci <- boot$total_summary
          names(total_ci)[names(total_ci) == "Boot_mean"] <- "Boot_mean_Total_effect"
          names(total_ci)[names(total_ci) == "Boot_SE"] <- "Boot_SE_Total_effect"
          names(total_ci)[names(total_ci) == "CI_low"] <- "Total_CI_low"
          names(total_ci)[names(total_ci) == "CI_high"] <- "Total_CI_high"
          names(total_ci)[names(total_ci) == "P_boot"] <- "P_boot_Total_effect"
          total_ci$Analysis <- analysis_label
          total_ci$Group <- group
          totals <- merge(
            totals, total_ci,
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
      lobo <- lobo_group(d, sp$equations, fit$edges)
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
    if (!is.null(cerror_edges) && nrow(cerror_edges)) {
      all_cerror_edges[[group]] <- cerror_edges
    }
    all_paths[[group]] <- paths
    all_totals[[group]] <- totals
    all_diag[[group]] <- diag
    all_fit[[group]] <- fit
  }

  combine_nonempty <- function(x) {
    x <- x[vapply(x, function(z) !is.null(z) && nrow(z) > 0, logical(1))]
    if (!length(x)) data.frame() else do.call(rbind, x)
  }

  edges <- combine_nonempty(all_edges)
  cerror_edges <- combine_nonempty(all_cerror_edges)
  paths <- combine_nonempty(all_paths)
  totals <- combine_nonempty(all_totals)
  diagnostics <- combine_nonempty(all_diag)
  lobo_stability <- combine_nonempty(all_lobo_stab)

  # Clean figure tables are constructed only after bootstrap and LOBO summaries
  # have been merged, so plotting code does not need to join diagnostic files.
  path_plot <- make_path_diagram_plot_table(
    edges = edges,
    correlated_error_edges = cerror_edges,
    lobo_stability = lobo_stability,
    alpha_display = alpha_display,
    figure_label = paste0(figure_label, "b")
  )
  effect_plot <- make_effect_plot_table(
    totals = totals,
    plot_source_order = plot_source_order,
    figure_label = paste0(figure_label, "c")
  )

  combined_dir <- file.path(out_dir, "combined")
  write_csv(edges, file.path(combined_dir, "01_all_fixed_path_coefficients.csv"))
  write_csv(
    cerror_edges,
    file.path(combined_dir, "01a_all_correlated_error_coefficients.csv")
  )
  write_csv(
    edges[edges$Display_significant, , drop = FALSE],
    file.path(combined_dir, paste0("02_", figure_label, "b_edges_significant_only.csv"))
  )
  write_csv(
    edges,
    file.path(combined_dir, paste0("03_", figure_label, "b_edges_all_fixed.csv"))
  )
  write_csv(
    path_plot,
    file.path(combined_dir, paste0("03a_", figure_label, "b_path_diagram_for_plot.csv"))
  )
  write_csv(paths, file.path(combined_dir, "04_all_unique_standardized_path_effects.csv"))
  write_csv(
    totals,
    file.path(combined_dir, paste0("05_", figure_label, "c_all_variable_total_effects.csv"))
  )
  write_csv(
    totals[totals$Source_variable %in% treatments, , drop = FALSE],
    file.path(combined_dir, paste0("05a_", figure_label, "c_treatment_total_effects.csv"))
  )
  # Primary Fig. c / Fig. S9c plotting input: every treatment and every
  # intermediate/environmental source is retained, including structural zeros.
  write_csv(
    effect_plot,
    file.path(combined_dir, paste0("06_", figure_label, "c_all_source_effects_for_plot.csv"))
  )
  write_csv(
    effect_plot[effect_plot$Source_type == "Treatment", , drop = FALSE],
    file.path(combined_dir, paste0("06a_", figure_label, "c_treatment_effects_for_plot.csv"))
  )
  write_csv(
    effect_plot[effect_plot$Source_variable == "warm", , drop = FALSE],
    file.path(combined_dir, paste0("06b_", figure_label, "c_warming_effects_for_plot.csv"))
  )
  write_csv(diagnostics, file.path(combined_dir, "07_component_model_diagnostics_R2.csv"))

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

  settings <- c(
    paste0("Analysis=", analysis_label),
    "Primary_method=mixed-effects piecewise SEM",
    paste0("piecewiseSEM_version=", as.character(utils::packageVersion("piecewiseSEM"))),
    paste0("lme4_version=", as.character(utils::packageVersion("lme4"))),
    "Path_selection=none within piecewiseSEM; MODEL_SPEC is fixed before fitting",
    paste0("Model_specification_basis=", model_spec_basis),
    "Component_model=response ~ fixed MODEL_SPEC parents + (1|block/WholePlotID)",
    "Random_structure=(1|block)+(1|block:WholePlotID)",
    "Block=random intercept matching the primary hierarchical analyses",
    "Block_WholePlotID=random intercept preserving clipped/unclipped subplot dependence",
    "Component_estimation=REML (required for piecewiseSEM F/Kenward-Roger tests)",
    "Model_call_environment=persistent per-component environment containing data and dat",
    "Coefficient=standardized path coefficient b*SD(x)/SD(y)",
    "Indirect_effect=product of standardized coefficients along a directed path",
    "Total_effect=direct effect + sum of all fixed indirect path effects",
    "Total_effect_sources=all non-target variables in each group MODEL_SPEC",
    "Variables_without_a_directed_path=retained with structural total effect 0",
    "Fig_c_for_plot=all treatments plus all intermediate/environmental variables, including structural zero effects",
    "Fig_b_for_plot=directed paths plus correlated errors in one clean table; Plot_include marks P<=Display_alpha",
    paste0("N_correlated_errors=", nrow(correlated_errors_by_group)),
    "Correlated_errors=passed to piecewiseSEM::psem as formula.cerror objects",
    "Correlated_errors_role=affect d-separation/Fisher's C; excluded from directed path and total-effect multiplication",
    paste0("Display_alpha=", alpha_display),
    "Display filtering does not alter the fitted model or total-effect calculation",
    paste0("Cluster_bootstrap=", run_bootstrap),
    paste0("N_boot=", if (run_bootstrap) n_boot else 0L),
    "Bootstrap_unit=whole plots resampled with replacement within each block",
    paste0("LOBO_sensitivity=", run_lobo),
    "LOBO_role=coefficient/sign robustness only; not model selection"
  )
  writeLines(settings, file.path(out_dir, "00_analysis_settings.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "00_sessionInfo.txt"))

  message(
    "\nFinished.\nClean path-diagram data: ",
    file.path(combined_dir, paste0("03a_", figure_label, "b_path_diagram_for_plot.csv")),
    "\nClean all-source effect data: ",
    file.path(combined_dir, paste0("06_", figure_label, "c_all_source_effects_for_plot.csv"))
  )

  invisible(list(
    data = dat,
    fits = all_fit,
    edges = edges,
    correlated_error_edges = cerror_edges,
    path_plot = path_plot,
    paths = paths,
    totals = totals,
    effect_plot = effect_plot,
    diagnostics = diagnostics
  ))
}
