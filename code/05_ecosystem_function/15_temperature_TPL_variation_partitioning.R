# Temperature and TPL variation partitioning, conditioning on block.
PROJECT_DIR <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
setwd(PROJECT_DIR)

# 1. USER SETTINGS ------------------------------------------------------------
ENV_FILE <-
  "./results/10_environment_SEM_inputs/all_environment_TPL_function_48plots.csv"
INFO_FILE <-
  "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"
OUTPUT_ROOT <- "./results/13_ecosystem_function/variation_partitioning"

ENV_PLOT_COLUMN <- NULL
INFO_PLOT_COLUMN <- NULL
INFO_BLOCK_COLUMN <- NULL
INFO_WHOLEPLOT_COLUMN <- NULL
INFO_CLIPPING_COLUMN <- NULL
INFO_WARMING_COLUMN <- NULL
INFO_PRECIPITATION_COLUMN <- NULL
RUN_TREATMENT_ADJUSTED <- FALSE  # retained only for the metadata reader

ENV_COLUMNS <- c("temperature_mean")
BETA_MAP <- data.frame(
  Group = rep(c("Bacteria", "Archaea", "Fungi", "Protists", "Nematodes"), 2),
  TPL = rep(c("tTTPL", "pTTPL"), each = 5),
  beta_variable = c("Bac_beta", "Arch_beta", "Fun_18S_beta", "Pro_beta", "Nem_beta",
                    "Bac_pTTPL_beta", "Arch_pTTPL_beta", "Fun_18S_pTTPL_beta",
                    "Pro_pTTPL_beta", "Nem_pTTPL_beta"),
  stringsAsFactors = FALSE)
FUNCTION_COLUMNS <- c(
  "GPP_mean", "GPP_cv", "NEE_mean", "NEE_cv", "ER_mean", "ER_cv",
  "Ra_mean", "Ra_cv", "Rh_mean", "Rh_cv", "Rs_mean", "Rs_cv")

RUN_MULTIVARIATE_RESPONSES <- TRUE
RUN_UNCONDITIONED_DESCRIPTIVE <- TRUE
MIN_RESIDUAL_DF <- 10L
N_PERM <- 9999L
SEED <- 20261014L
RESPONSE_TRANSFORM <- "none"  # alternative: "log1p_nonnegative"
INSTALL_MISSING_PACKAGES <- FALSE
WRITE_MODEL_RDS <- TRUE
MAKE_PLOTS <- TRUE

# Primary multiplicity adjustment is BH within response family AND hypothesis,
# across all ten exponents. For the default responses this is 120 tests for
# TPL_unique on single functions and 30 for TPL_unique on multivariate profiles.
# Temperature_unique and joint tests form separate, equally defined families.
# Families do not reset for individual groups, exponents or individual functions.
# q_BH_global is an additional adjustment across all 450 permutation tests.

for (nm in c("ENV_FILE", "INFO_FILE", "OUTPUT_ROOT")) {
  val <- Sys.getenv(paste0("TPL_", nm), unset = "")
  if (nzchar(val)) assign(nm, val)
}
if (nzchar(Sys.getenv("TPL_N_PERM"))) N_PERM <- as.integer(Sys.getenv("TPL_N_PERM"))

main <- function() {
  if (!RESPONSE_TRANSFORM %in% c("none", "log1p_nonnegative"))
    stop("RESPONSE_TRANSFORM must be none or log1p_nonnegative.")
  if (!is.finite(N_PERM) || N_PERM < 19L) stop("N_PERM must be >= 19.")
  if (!identical(ENV_COLUMNS, "temperature_mean"))
    stop("This focused analysis requires ENV_COLUMNS <- c('temperature_mean').")
  if (anyDuplicated(BETA_MAP$beta_variable)) stop("Duplicated entries in BETA_MAP.")
  pkgs <- c("vegan", "permute")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && INSTALL_MISSING_PACKAGES)
    install.packages(missing, repos = "https://cloud.r-project.org")
  if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
    stop("Install vegan and permute, then rerun.")
  paths <- c(ENV_FILE, INFO_FILE)
  if (any(!file.exists(paths)))
    stop("Missing input file(s):\n", paste(paths[!file.exists(paths)], collapse = "\n"),
         "\nWorking directory: ", getwd())
  set.seed(SEED)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(OUTPUT_ROOT, paste0("run_", stamp))
  if (dir.exists(out)) out <- paste0(out, "_", Sys.getpid())
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (d in c("input_checks", "models", "figures")) dir.create(file.path(out, d))
  logcon <- file(file.path(out, "run_log.txt"), "wt")
  sink(logcon, split = TRUE)
  on.exit({sink(); close(logcon)}, add = TRUE)
  cat("Output:", normalizePath(out), "\n")
  write_csv <- function(x, p) utils::write.csv(x, file.path(out, p),
      row.names = FALSE, na = "", fileEncoding = "UTF-8")
  norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  norm_id <- function(x) {
    x <- toupper(trimws(as.character(x)))
    # N01_2010 / S1-2010 -> N01 / S01, for annual treatment metadata only.
    m <- regexec("^([NS])0*([0-9]{1,2})($|[_. -])", x)
    hits <- regmatches(x, m)
    ok <- lengths(hits) > 0L
    x[ok] <- vapply(hits[ok], function(z) paste0(z[2], sprintf("%02d", as.integer(z[3]))), "")
    x
  }
  read_tab <- function(path) {
    d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA", "NaN", "NULL"), fileEncoding = "UTF-8-BOM")
    blank <- !nzchar(names(d))
    names(d)[blank] <- paste0("Unnamed_column_", which(blank))
    if (!nrow(d) || !ncol(d)) stop("Empty table: ", path)
    if (anyDuplicated(names(d))) stop("Duplicated column names: ", path)
    d
  }
  choose_id <- function(d, specified, ref = NULL) {
    if (!is.null(specified)) {
      if (!specified %in% names(d)) stop("ID column not found: ", specified)
      return(specified)
    }
    if (is.null(ref)) {
      preferred <- which(norm_name(names(d)) %in% c("plot", "plotid", "subplot", "subplotid", "sampleid"))
      i <- if (length(preferred)) preferred[1] else 1L
    } else {
      scores <- vapply(d, function(z) length(intersect(unique(norm_id(z)), ref)), integer(1))
      i <- which.max(scores)
      if (scores[i] != length(ref))
        stop("Cannot identify a column matching all ENV plot IDs. Set the corresponding *_PLOT_COLUMN.")
    }
    names(d)[i]
  }
  keyed <- function(d, key, label) {
    ids <- norm_id(d[[key]])
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop(label, ": blank or duplicated plot IDs.")
    d[[key]] <- NULL
    rownames(d) <- ids
    d
  }
  check_numeric <- function(d, label) {
    if (!ncol(d)) stop("No numeric predictors in ", label)
    for (j in seq_along(d)) {
      old <- d[[j]]
      num <- suppressWarnings(as.numeric(as.character(old)))
      if (any(!is.na(old) & is.na(num)))
        stop(label, ": nonnumeric values in ", names(d)[j], ". Only the plot-ID column may be text.")
      num[!is.finite(num)] <- NA_real_
      d[[j]] <- num
    }
    d
  }

  env0 <- read_tab(ENV_FILE)
  env <- keyed(env0, choose_id(env0, ENV_PLOT_COLUMN), "ENV_FILE")
  ids <- rownames(env)
  required <- unique(c(ENV_COLUMNS, BETA_MAP$beta_variable, FUNCTION_COLUMNS))
  if (any(!required %in% names(env)))
    stop("ENV missing columns: ", paste(setdiff(required, names(env)), collapse = ", "))
  env <- check_numeric(env[required], "Selected ENV columns")
  E0 <- env[ENV_COLUMNS]
  B0 <- env[BETA_MAP$beta_variable]
  Y0 <- env[FUNCTION_COLUMNS]
  info <- NULL
  if (file.exists(INFO_FILE)) {
    raw <- read_tab(INFO_FILE)
    ikey <- choose_id(raw, INFO_PLOT_COLUMN, ids)
    pids <- norm_id(raw[[ikey]])
    pick <- function(specified, aliases, required = TRUE) {
      if (!is.null(specified)) {
        if (!specified %in% names(raw)) stop("INFO column not found: ", specified)
        return(specified)
      }
      hits <- names(raw)[norm_name(names(raw)) %in% aliases]
      if (length(hits)) return(hits[1])
      if (required) stop("Cannot find INFO column for ", aliases[1],
                         ". Set the corresponding INFO_*_COLUMN explicitly.")
      NULL
    }
    maps <- list(Block = pick(INFO_BLOCK_COLUMN, c("block", "blockid")),
                 WholePlotID = pick(INFO_WHOLEPLOT_COLUMN, c("wholeplotid", "wholeplot")),
                 Clipping = pick(INFO_CLIPPING_COLUMN, c("clipping", "clip", "clipped", "mowing")),
                 Warming = pick(INFO_WARMING_COLUMN, c("warming", "warm", "temperaturetreatment"), RUN_TREATMENT_ADJUSTED),
                 Precipitation = pick(INFO_PRECIPITATION_COLUMN,
                                      c("precipitation", "precip", "ppt", "precipitationtreatment"), RUN_TREATMENT_ADJUSTED))
    info <- data.frame(Plot = ids, stringsAsFactors = FALSE)
    for (nm in names(maps)) {
      col <- maps[[nm]]
      if (is.null(col)) next
      values <- vapply(ids, function(id) {
        z <- unique(trimws(as.character(raw[pids == id, col])))
        if (length(z) != 1 || anyNA(z) || !nzchar(z))
          stop("Missing/inconsistent INFO ", nm, " within ", id,
               ". Annual duplicates must have identical treatment metadata.")
        z
      }, "")
      info[[nm]] <- factor(values)
    }
    info$Cluster <- interaction(info$Block, info$WholePlotID, drop = TRUE)
    rownames(info) <- ids
    map_table <- data.frame(role = c("Plot", names(maps)),
                            source_column = c(ikey, vapply(maps, function(z) if (is.null(z)) "not_used" else z, "")))
    write_csv(map_table, "input_checks/metadata_column_map.csv")
  } else cat("INFO absent: descriptive mode. No inferential P-values.\n")


  # Common complete cases keep all comparisons on the same subplot sample set.
  # If one subplot is missing, remove its entire paired whole plot.
  finite <- stats::complete.cases(E0, B0, Y0)
  reasons <- ifelse(finite, "retained", "missing_selected_value")
  bad <- unique(info$Cluster[!finite])
  paired_bad <- info$Cluster %in% bad & finite
  finite[paired_bad] <- FALSE
  reasons[paired_bad] <- "paired_wholeplot_removed"
  write_csv(data.frame(Plot = ids, included = finite, reason = reasons),
            "input_checks/sample_inclusion.csv")
  ids <- ids[finite]
  info <- droplevels(info[ids, , drop = FALSE])
  info <- info[order(info$Block, info$Cluster, info$Clipping), , drop = FALSE]
  ids <- rownames(info)
  E0 <- E0[ids, , drop = FALSE]
  B0 <- B0[ids, , drop = FALSE]
  Y0 <- Y0[ids, , drop = FALSE]
  n <- length(ids)
  if (n < 20L) stop("Fewer than 20 complete subplot trajectories remain.")
  if (nlevels(info$Block) < 2L) stop("At least two blocks are required.")
  vary <- vapply(cbind(E0, B0, Y0), stats::sd, numeric(1))
  if (any(!is.finite(vary) | vary <= 1e-12))
    stop("Non-varying or invalid selected columns: ",
         paste(names(vary)[!is.finite(vary) | vary <= 1e-12], collapse = ", "))
  write_csv(data.frame(Plot = ids, E0, B0, Y0, check.names = FALSE),
            "input_checks/aligned_analysis_data.csv")
  write_csv(info, "input_checks/aligned_metadata.csv")
  write_csv(BETA_MAP, "input_checks/beta_variable_map.csv")
  write_csv(data.frame(input_file = paths, md5 = unname(tools::md5sum(paths))),
            "input_checks/input_md5.csv")
  cat("Complete subplots:", n, "; blocks:", nlevels(info$Block),
      "; whole plots:", nlevels(info$Cluster), "\n")
  E <- as.matrix(scale(E0))
  B <- as.matrix(scale(B0))
  response_transforms <- data.frame(variable = names(Y0), transform = "none")
  if (RESPONSE_TRANSFORM == "log1p_nonnegative") for (j in seq_along(Y0)) {
    if (all(Y0[[j]] >= 0)) {
      Y0[[j]] <- log1p(Y0[[j]])
      response_transforms$transform[j] <- "log1p"
    } else response_transforms$transform[j] <- "none_signed_values"
  }
  write_csv(response_transforms, "input_checks/response_transformations.csv")
  Y <- as.matrix(scale(Y0))
  responses <- setNames(lapply(FUNCTION_COLUMNS, function(z) z), FUNCTION_COLUMNS)
  if (RUN_MULTIVARIATE_RESPONSES) {
    responses$ALL_12_FUNCTIONS <- FUNCTION_COLUMNS
    responses$MEAN_6_FUNCTIONS <- FUNCTION_COLUMNS[grepl("_mean$", FUNCTION_COLUMNS)]
    responses$CV_6_FUNCTIONS <- FUNCTION_COLUMNS[grepl("_cv$", FUNCTION_COLUMNS)]
  }
  # Block_adjusted is the primary analysis. Unconditioned fractions are descriptive.
  Cblock <- stats::model.matrix(~ Block, info)[, -1, drop = FALSE]
  conditions <- list(Block_adjusted = Cblock)
  if (RUN_UNCONDITIONED_DESCRIPTIVE) conditions["unconditioned"] <- list(NULL)
  permutation_matrix <- NULL
  if (!is.null(info)) {
    signatures <- lapply(split(as.character(info$Clipping), info$Cluster), sort)
    if (length(unique(vapply(signatures, paste, "", collapse = ":"))) != 1L ||
        any(vapply(signatures, anyDuplicated, integer(1)) > 0L))
      stop("Whole plots must have the same unique Clipping positions for paired permutations.")
    ctrl <- permute::how(nperm = N_PERM, blocks = info$Block,
                         plots = permute::Plots(strata = info$Cluster, type = "free"),
                         within = permute::Within(type = "none"))
    permutation_matrix <- permute::shuffleSet(n, nset = N_PERM, control = ctrl)
    # Assert that permutations never cross blocks or clipping positions.
    for (i in seq_len(nrow(permutation_matrix))) {
      p <- permutation_matrix[i, ]
      stopifnot(all(info$Block[p] == info$Block), all(info$Clipping[p] == info$Clipping))
      mapped_clusters <- tapply(as.character(info$Cluster[p]), info$Cluster,
                                function(z) length(unique(z)))
      stopifnot(all(mapped_clusters == 1L))
    }
    saveRDS(permutation_matrix, file.path(out, "input_checks/permutation_matrix.rds"))
    cat("Permutations:", nrow(permutation_matrix),
        "; whole plots shuffled within Block, paired Clipping order preserved.\n")
  }

  rank_centered <- function(x) {
    if (is.null(x) || !ncol(x)) return(0L)
    qr(scale(x, center = TRUE, scale = FALSE), tol = 1e-8)$rank
  }
  fit <- function(y, x, z = NULL) {
    if (is.null(z) || !ncol(z)) vegan::rda(y, x) else vegan::rda(y, x, z)
  }
  rsq <- function(y, x) {
    if (is.null(x) || !ncol(x)) return(c(raw = 0, adjusted = 0))
    rr <- vegan::RsquareAdj(fit(y, x))
    c(raw = rr$r.squared, adjusted = rr$adj.r.squared)
  }
  permutation_test <- function(model) {
    if (is.null(model$CCA) || model$CCA$rank == 0L)
      return(data.frame(F = NA_real_, df = 0L, p = NA_real_, permutations = 0L,
                        test_status = "no_unique_predictor_rank"))
    aa <- stats::anova(model, permutations = permutation_matrix,
                       model = "reduced", parallel = 1L)
    data.frame(F = aa$F[1], df = aa$Df[1], p = aa$`Pr(>F)`[1],
               permutations = nrow(permutation_matrix), test_status = "tested")
  }
  results <- list(); fractions <- list(); tests <- list(); diagnostics <- list()
  nr <- 0L; nt <- 0L; nd <- 0L
  comps <- c("Temperature_unique", "TPL_unique", "Shared", "Unexplained")
  for (i in seq_len(nrow(BETA_MAP))) for (baseline in names(conditions)) {
    spec <- BETA_MAP[i, , drop = FALSE]
    xB <- B[, spec$beta_variable, drop = FALSE]
    C <- conditions[[baseline]]
    rFull <- rank_centered(cbind(C, E, xB))
    rdf <- n - 1L - rFull
    strategy <- paste(spec$Group, spec$TPL, sep = "_")
    for (response in names(responses)) {
      cat("Analysis:", strategy, "/", baseline, "/", response, "\n")
      nd <- nd + 1L
      ident <- data.frame(strategy = strategy, Group = spec$Group, TPL = spec$TPL,
                          beta_variable = spec$beta_variable, baseline = baseline,
                          response = response, stringsAsFactors = FALSE)
      diagnostics[[nd]] <- cbind(ident, n = n, joint_predictor_rank = rFull,
          condition_rank = rank_centered(C), residual_df = rdf,
          status = "pending", message = "")
      if (rdf < MIN_RESIDUAL_DF) {
        diagnostics[[nd]]$status <- "skipped_insufficient_df"
        diagnostics[[nd]]$message <- paste("Residual df", rdf, "<", MIN_RESIDUAL_DF)
        next
      }
      y <- Y[, responses[[response]], drop = FALSE]
      tryCatch({
        RC <- rsq(y, C)
        REC <- rsq(y, cbind(C, E))
        RBC <- rsq(y, cbind(C, xB))
        REBC <- rsq(y, cbind(C, E, xB))
        partition <- function(j) c(
          Temperature_unique = REBC[j] - RBC[j],
          TPL_unique = REBC[j] - REC[j],
          Shared = REC[j] + RBC[j] - REBC[j] - RC[j],
          Unexplained = 1 - REBC[j])
        rawf <- unname(partition("raw"))
        adjf <- unname(partition("adjusted"))
        names(rawf) <- names(adjf) <- comps
        remaining <- unname(1 - RC["adjusted"])
        if (!is.finite(remaining) || remaining <= 1e-8)
          stop("No remaining adjusted variation after block covariates.")
        stopifnot(abs(sum(adjf) + RC["adjusted"] - 1) < 1e-7)
        if (is.null(C)) {
          vp <- vegan::varpart(y, E, xB)
          stopifnot(max(abs(vp$part$indfract$Adj.R.squared[1:4] - adjf)) < 1e-7)
        } else vp <- NULL
        negative <- any(adjf[1:3] < -1e-8)
        nr <- nr + 1L
        results[[nr]] <- cbind(ident,
          response_kind = if (ncol(y) == 1L) "single" else "multivariate",
          n = n, response_columns = ncol(y), env_dimensions = ncol(E),
          beta_dimensions = ncol(xB), residual_df = rdf,
          covariate_R2_raw = unname(RC["raw"]),
          covariate_R2_adjusted = unname(RC["adjusted"]),
          env_plus_covariates_R2_adjusted = unname(REC["adjusted"]),
          beta_plus_covariates_R2_adjusted = unname(RBC["adjusted"]),
          full_R2_raw = unname(REBC["raw"]), full_R2_adjusted = unname(REBC["adjusted"]),
          env_total_remaining = unname((REC["adjusted"] - RC["adjusted"]) / remaining),
          TPL_total_remaining = unname((RBC["adjusted"] - RC["adjusted"]) / remaining),
          joint_explained_remaining = unname((REBC["adjusted"] - RC["adjusted"]) / remaining),
          negative_fraction_present = negative)
        fractions[[nr]] <- cbind(ident[rep(1L, 4L), , drop = FALSE],
          component = comps, fraction_total_raw = rawf,
          fraction_total_adjusted = adjf,
          fraction_remaining_adjusted = adjf / remaining,
          percent_remaining_adjusted = 100 * adjf / remaining,
          denominator_remaining_adjusted = remaining)
        mods <- list(joint = fit(y, cbind(E, xB), C),
          Temperature_unique = fit(y, E, cbind(C, xB)),
          TPL_unique = fit(y, xB, cbind(C, E)))
        if (baseline == "Block_adjusted") for (h in names(mods)) {
          nt <- nt + 1L
          tests[[nt]] <- cbind(ident, hypothesis = h, permutation_test(mods[[h]]))
        }
        if (WRITE_MODEL_RDS) saveRDS(list(models = mods, varpart = vp,
            fractions = fractions[[nr]], covariates = C,
            temperature = E, beta = xB, response = y),
            file.path(out, "models", paste0(strategy, "__", baseline, "__", response, ".rds")))
        diagnostics[[nd]]$status <- "completed"
        diagnostics[[nd]]$message <- if (negative) "Signed negative adjusted fractions retained." else ""
      }, error = function(e) {
        diagnostics[[nd]]$status <<- "error"
        diagnostics[[nd]]$message <<- conditionMessage(e)
        cat("ERROR:", conditionMessage(e), "\n")
      })
    }
  }
  status <- do.call(rbind, diagnostics)
  write_csv(status, "analysis_status.csv")
  if (!length(results)) stop("No analysis completed; inspect analysis_status.csv.")
  res <- do.call(rbind, results)
  fr <- do.call(rbind, fractions)
  write_csv(res, "model_summary.csv")
  write_csv(fr, "variation_fractions.csv")
  # Two-predictor hierarchical partitioning: average each predictor's contribution
  # over both possible entry orders, retaining signed adjusted fractions.
  hp <- fr[fr$component == "TPL_unique", , drop = FALSE]
  shared <- fr[fr$component == "Shared", , drop = FALSE]
  temp <- fr[fr$component == "Temperature_unique", , drop = FALSE]
  hp$TPL_independent_adjusted <- hp$fraction_remaining_adjusted + shared$fraction_remaining_adjusted / 2
  hp$Temperature_independent_adjusted <- temp$fraction_remaining_adjusted + shared$fraction_remaining_adjusted / 2
  write_csv(hp[, c("Group", "TPL", "beta_variable", "baseline", "response",
                   "TPL_independent_adjusted", "Temperature_independent_adjusted")],
            "hierarchical_partitioning.csv")
  ts <- if (length(tests)) do.call(rbind, tests) else data.frame()
  if (nrow(ts)) {
    ts$response_family <- ifelse(ts$response %in% FUNCTION_COLUMNS,
                                 "single_functions", "multivariate_profiles")
    ts$correction_family <- paste(ts$baseline, ts$response_family, ts$hypothesis, sep = "__")
    ts$q_BH_family <- NA_real_
    ts$family_size <- NA_integer_
    for (ix in split(seq_len(nrow(ts)), ts$correction_family)) {
      ts$q_BH_family[ix] <- stats::p.adjust(ts$p[ix], method = "BH")
      ts$family_size[ix] <- length(ix)
    }
    ts$q_BH_global <- stats::p.adjust(ts$p, method = "BH")
    write_csv(ts, "permutation_tests.csv")
    fam <- unique(ts[c("correction_family", "response_family", "hypothesis", "family_size")])
    write_csv(fam, "multiple_testing_families.csv")
  } else write_csv(data.frame(note = "No permutation tests completed."), "permutation_tests.csv")

  # One row per group x TPL type x response, for the PRIMARY block-adjusted model.
  b16 <- res[res$baseline == "Block_adjusted", , drop = FALSE]
  keys <- c("strategy", "Group", "TPL", "beta_variable", "baseline", "response")
  for (comp in comps) {
    dd <- fr[fr$baseline == "Block_adjusted" & fr$component == comp,
             c(keys, "fraction_total_raw", "fraction_total_adjusted", "percent_remaining_adjusted")]
    names(dd)[!names(dd) %in% keys] <- paste0(comp, c("_raw", "_adjusted", "_percent_remaining"))
    b16 <- merge(b16, dd, by = keys, all.x = TRUE, sort = FALSE)
  }
  if (nrow(ts)) for (h in c("TPL_unique", "Temperature_unique", "joint")) {
    dd <- ts[ts$hypothesis == h, c(keys, "F", "df", "p", "q_BH_family", "family_size", "q_BH_global")]
    names(dd)[!names(dd) %in% keys] <- paste0(h, "_", c("F", "df", "p", "q_BH_family", "family_size", "q_BH_global"))
    b16 <- merge(b16, dd, by = keys, all.x = TRUE, sort = FALSE)
  }
  b16 <- b16[order(match(b16$beta_variable, BETA_MAP$beta_variable),
                    match(b16$response, names(responses))), , drop = FALSE]
  write_csv(b16, "B16_single_beta_summary.csv")

  if (MAKE_PLOTS) for (strategy in unique(fr$strategy)) {
    d <- fr[fr$strategy == strategy & fr$baseline == "Block_adjusted", , drop = FALSE]
    if (!nrow(d)) next
    grDevices::pdf(file.path(out, "figures", paste0(strategy, "__fractions.pdf")),
                   width = 11, height = 8.5, onefile = TRUE)
    op <- par(mfrow = c(3, 3), mar = c(4, 8, 3, 1), oma = c(2, 0, 2, 0))
    cols <- c("#4C78A8", "#F58518", "#54A24B", "#B8B8B8")
    for (response in names(responses)) {
      dd <- d[d$response == response, ]
      if (!nrow(dd)) next
      v <- dd$percent_remaining_adjusted[match(comps, dd$component)]
      lim <- range(c(v, 0)); pad <- max(diff(lim), 1) * 0.20
      pos <- barplot(rev(v), horiz = TRUE, col = rev(cols), border = NA,
          names.arg = rev(c("Temperature only", "TPL beta only", "Shared", "Unexplained")),
          las = 1, cex.names = 0.75, xlim = c(lim[1] - pad, lim[2] + pad),
          main = paste(strategy, response, sep = ": "), cex.main = 0.8,
          xlab = "% of variation remaining after block adjustment")
      abline(v = 0, col = "#555555")
      text(rev(v), pos, sprintf("%.1f", rev(v)),
           pos = ifelse(rev(v) >= 0, 4, 2), cex = 0.75)
    }
    par(op); grDevices::dev.off()
  }
  settings <- list(input_paths = c(ENV_FILE = ENV_FILE, INFO_FILE = INFO_FILE),
      environment_columns = ENV_COLUMNS, beta_map = BETA_MAP,
      function_columns = FUNCTION_COLUMNS, responses = responses,
      seed = SEED, permutations = nrow(permutation_matrix),
      min_residual_df = MIN_RESIDUAL_DF, transform = RESPONSE_TRANSFORM,
      plot_ids = ids, primary_baseline = "Block_adjusted",
      design = "Whole plots permuted within block; clipping positions preserved",
      multiplicity = "BH across 10 exponents and all responses within each response family and hypothesis")
  saveRDS(settings, file.path(out, "analysis_settings.rds"))
  capture.output(str(settings), file = file.path(out, "analysis_settings.txt"))
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  writeLines(c(
    "READ FIRST: mean soil temperature vs one TPL exponent",
    "Ten separate two-set partitions: five groups x tTTPL/pTTPL; no PCA or subgroup-file inputs.",
    "B16_single_beta_summary.csv: PRIMARY Block_adjusted fractions and tests joined by analysis.",
    "variation_fractions.csv: original and adjusted fractions; negative values are retained.",
    "For Block_adjusted, fraction_total_adjusted is on the original response-variation scale.",
    "The four fraction_total_adjusted values plus covariate_R2_adjusted sum to one.",
    "percent_remaining_adjusted = 100 * fraction_total_adjusted / (1 - covariate_R2_adjusted).",
    "These four normalized percentages sum to 100. Negative block-adjusted R2 can make their denominator exceed one.",
    "Unconditioned partitions are descriptive and have no permutation P-values.",
    "TPL_unique tests ONE exponent after mean soil temperature and block, not after other groups/exponents.",
    "Temperature_unique tests mean soil temperature after that exponent and block.",
    "Shared fractions are overlaps, not evidence for causality, mediation or interaction; they have no direct P-value.",
    "q_BH_family: within each hypothesis, adjust across all ten exponents and all responses of the same family.",
    "Default family sizes: 120 single-function tests or 30 multivariate-profile tests for EACH hypothesis.",
    "q_BH_global: sensitivity adjustment across all 450 permutation tests (default settings).",
    "No per-group or per-function resetting of the primary BH families is performed.",
    "Beta t/p are tested separately; these partitions do not isolate beta_t from beta_p.",
    "Permutation inference assumes paired residual profiles are exchangeable among whole plots within blocks.",
    "Multivariate responses are standardized function profiles, not multifunctionality indices.",
    "NEE CV can be unstable when its temporal mean is close to zero; check the source time series.",
    "Check analysis_status.csv before reporting any results."
  ), file.path(out, "READ_FIRST.txt"))
  cat("\nCompleted:", sum(status$status == "completed"),
      "; skipped:", sum(grepl("^skipped", status$status)),
      "; errors:", sum(status$status == "error"), "\n")
  if (any(status$status == "error")) warning("Some analyses failed; inspect analysis_status.csv.")
  cat("Results saved to:", normalizePath(out), "\n")
  invisible(list(output_directory = out, model_summary = res,
                  fractions = fr, permutation_tests = ts, B16_summary = b16, status = status))
}

TPL_VARIATION_RESULTS <- main()
