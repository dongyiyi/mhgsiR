#' Diagnose a Symbiotic Mechanism from a Fingerprint
#'
#' @description
#' Compare an observed macro‑fingerprint to a library of theoretical archetypes
#' and rank the best matches.
#'
#' @details
#' This is the core inferential function of the mhgsiR package. It takes a
#' 5D fingerprint (e.g., from `extract_fingerprint()`) and computes a
#' similarity score against each reference fingerprint. The choice of method and
#' rescaling can influence results, and users are encouraged to test for concordance
#' across different methods.
#'
#' ### Choosing a Similarity Method:
#'
#' * **`cosine` (Default):** Measures the orientation (angle) of fingerprint vectors.
#'   - **Use when:** The relative pattern and direction of effects are most important. Robust for low-dimensional data (≤3 metrics).
#'   - **Use with caution when:** The absolute magnitude of differences is biologically crucial.
#'
#' * **`correlation` (Pearson):** Measures the linear relationship of the fingerprint profiles (shape).
#'   - **Use when:** You care about the "shape" of the profile, and metrics have been rescaled to a comparable range. More robust with ≥4 metrics.
#'   - **Use with caution when:** The absolute values or magnitudes are biologically meaningful, as this information is weakened by mean-centering.
#'
#' * **`euclidean`:** Measures the absolute distance between fingerprints in multi-dimensional space.
#'   - **Use when:** The goal is to find the closest match based on absolute numerical values. Requires metrics to be on a comparable scale (i.e., `allow_rescale = TRUE` is highly recommended).
#'   - **Use with caution when:** Metrics have vastly different scales and rescaling is turned off, as the metric with the largest scale will dominate the distance.
#'
#' ### Choosing to Rescale (`allow_rescale`):
#'
#' * **`FALSE` (Default):** Compares raw metric values. Preserves all information about the magnitude of effects.
#'   - **Recommended for:** Analyzing data from a single, consistent experiment where the scales of the metrics are directly meaningful and comparable.
#'
#' * **`TRUE`:** Applies a z-score transformation to each metric across all fingerprints being compared.
#'   - **Recommended for:** Comparing fingerprints from different experiments, studies, or systems where the original scales may not be directly comparable. This focuses the comparison on "shape" rather than "magnitude".
#'
#' @param observed_fingerprint A 1-row tibble from `extract_fingerprint()`.
#' @param reference_library A data.frame of reference fingerprints. Defaults to `mhgsiR::reference_fingerprints`.
#' @param method The similarity metric. Defaults to "cosine".
#' @param allow_rescale Logical. If FALSE (default), compares raw values. If TRUE, applies z-score scaling.
#' @param na_strategy "drop_metric" (default) drops metrics with NA/NaN/Inf in either side;
#'   "error" stops on any NA/NaN/Inf.
#' @param weights Optional numeric vector (only for "weighted_correlation").
#'
#' @return A tibble of class `mhgsi_diagnosis` with columns:
#'   Mechanism, similarity_score, method, n_metrics, rank, tie_rank.
#'   Attributes: used_metrics, observed_vector, metric_names_used.
#' @export
diagnose_mechanism <- function(
    observed_fingerprint,
    reference_library = mhgsiR::reference_fingerprints,
    method = c("cosine", "correlation", "euclidean", "weighted_correlation"),
    allow_rescale = FALSE,
    na_strategy = c("drop_metric", "error"),
    weights = NULL
) {
  method <- match.arg(method)
  na_strategy <- match.arg(na_strategy)
  metric_names <- .mhgsi_metric_names

  # ---- compat: accept legacy mechanism column names ----
  if (!"Mechanism" %in% names(reference_library)) {
    alt <- c("mechanism_observed","mechanism","mechanism_name")
    hit <- alt[alt %in% names(reference_library)][1]
    if (!is.na(hit)) {
      reference_library <- dplyr::rename(reference_library, Mechanism = !!rlang::sym(hit))
    } else {
      stop("`reference_library` must contain a `Mechanism` column.", call. = FALSE)
    }
  }

  # ---- validation ----
  if (!is.data.frame(observed_fingerprint) || nrow(observed_fingerprint) != 1L)
    stop("`observed_fingerprint` must be a data.frame/tibble with exactly one row.", call. = FALSE)
  if (!is.data.frame(reference_library) || nrow(reference_library) < 1L)
    stop("`reference_library` must be a non-empty data.frame/tibble.", call. = FALSE)
  if (!all(metric_names %in% names(observed_fingerprint)) ||
      !all(metric_names %in% names(reference_library)))
    stop("Both inputs must contain all 5 core metric columns.", call. = FALSE)

  # ---- extract & coerce to numeric ----
  observed_vec    <- as.numeric(as.vector(unlist(observed_fingerprint[1, metric_names])))
  reference_matrix <- apply(reference_library[, metric_names, drop = FALSE], 2, as.numeric)

  # ---- NA/NaN/Inf filtering (moved earlier; governs both sides) ----
  is_finite_obs     <- is.finite(observed_vec)
  is_finite_ref_col <- apply(is.finite(reference_matrix), 2, all)

  if (na_strategy == "drop_metric") {
    keep_cols <- is_finite_obs & is_finite_ref_col
    if (!any(keep_cols))
      stop("All metrics are NA/NaN/Inf after filtering; cannot compute similarity.", call. = FALSE)

    observed_vec_clean     <- observed_vec[keep_cols]
    reference_matrix_clean <- reference_matrix[, keep_cols, drop = FALSE]
    used_metrics           <- metric_names[keep_cols]
    dropped_na_or_inf      <- metric_names[!keep_cols]
  } else { # na_strategy == "error"
    if (!all(is_finite_obs))
      stop("Observed fingerprint contains non-finite values.", call. = FALSE)
    if (!all(is_finite_ref_col))
      stop("Reference library contains non-finite values in some metric columns.", call. = FALSE)
    observed_vec_clean     <- observed_vec
    reference_matrix_clean <- reference_matrix
    used_metrics           <- metric_names
    dropped_na_or_inf      <- character(0)
  }

  # ---- drop zero-variance metrics (before scaling) ----
  zero_var <- apply(rbind(observed_vec_clean, reference_matrix_clean), 2, function(v) stats::sd(v) == 0)
  if (any(zero_var)) {
    observed_vec_clean     <- observed_vec_clean[!zero_var]
    reference_matrix_clean <- reference_matrix_clean[, !zero_var, drop = FALSE]
    used_metrics           <- used_metrics[!zero_var]
  }
  dropped_zero <- setdiff(metric_names, used_metrics)

  # ---- message which metrics were dropped ----
  if (length(dropped_na_or_inf) || length(dropped_zero)) {
    msg <- c()
    if (length(dropped_na_or_inf)) msg <- c(msg, paste0("NA/NaN/Inf: ", paste(dropped_na_or_inf, collapse = ", ")))
    if (length(dropped_zero))      msg <- c(msg, paste0("zero-variance: ", paste(dropped_zero, collapse = ", ")))
    message("Dropped metrics -> ", paste(msg, collapse = " | "))
  }
  if (length(used_metrics) == 0L)
    stop("No metrics left after filtering.", call. = FALSE)

  # ---- optional z-score across obs+ref per metric ----
  if (isTRUE(allow_rescale)) {
    cm <- rbind(observed_vec_clean, reference_matrix_clean)
    sm <- scale(cm)
    sm[is.nan(sm)] <- 0
    observed_vec_clean     <- sm[1, ]
    reference_matrix_clean <- sm[-1, , drop = FALSE]
  }

  # ---- weights (for weighted_correlation) ----
  if (identical(method, "weighted_correlation")) {
    if (is.null(weights))
      stop("`weights` must be provided for method = 'weighted_correlation'.", call. = FALSE)
    if (length(weights) != length(used_metrics))
      stop("`weights` length must equal number of used metrics (", length(used_metrics), ").", call. = FALSE)
    if (any(weights < 0) || sum(weights) == 0)
      stop("`weights` must be non-negative and sum to > 0.", call. = FALSE)
    w <- as.numeric(weights) / sum(weights)
  }

  # ---- similarity functions (guards included) ----
  pearson_cor <- function(a, b) {
    if (stats::sd(a) == 0 || stats::sd(b) == 0) return(0)
    val <- suppressWarnings(stats::cor(a, b, method = "pearson"))
    ifelse(is.na(val), 0, val)
  }
  cosine_sim <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (d == 0) return(0)
    sum(a * b) / d
  }
  euclid_sim <- function(a, b) -sqrt(sum((a - b)^2))
  wcor <- function(a, b, w) {
    wa <- sum(w * a); wb <- sum(w * b)
    num <- sum(w * (a - wa) * (b - wb))
    den <- sqrt(sum(w * (a - wa)^2)) * sqrt(sum(w * (b - wb)^2))
    if (den == 0) return(0)
    num / den
  }

  sim_fun <- switch(method,
                    correlation          = pearson_cor,
                    cosine               = cosine_sim,
                    euclidean            = euclid_sim,
                    weighted_correlation = function(a, b) wcor(a, b, w)
  )

  scores <- apply(reference_matrix_clean, 1, function(r) sim_fun(observed_vec_clean, r))

  # ---- output ----
  if (anyDuplicated(reference_library$Mechanism))
    message("Note: `reference_library$Mechanism` has duplicates; returning per-row scores.")

  result <- tibble::tibble(
    Mechanism = reference_library$Mechanism,
    similarity_score = as.numeric(scores),
    method = method,
    n_metrics = length(used_metrics)
  ) |>
    dplyr::arrange(dplyr::desc(similarity_score)) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      tie_rank = dplyr::min_rank(dplyr::desc(similarity_score))
    )

  attr(result, "used_metrics") <- used_metrics
  attr(result, "observed_vector") <- as.numeric(observed_vec_clean)
  attr(result, "metric_names_used") <- used_metrics
  class(result) <- c("mhgsi_diagnosis", class(result))
  result
}


#' @export
#' @rdname diagnose_mechanism
plot.mhgsi_diagnosis <- function(x, ...) {
  ggplot2::ggplot(
    x, ggplot2::aes(x = similarity_score, y = stats::reorder(Mechanism, similarity_score))
  ) +
    ggplot2::geom_col(ggplot2::aes(fill = rank == 1L), alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", similarity_score)),
      hjust = -0.15, size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey70"), guide = "none") +
    ggplot2::labs(
      title = "Mechanism Diagnosis",
      subtitle = sprintf("Method: %s • Metrics used: %d", x$method[1], x$n_metrics[1]),
      x = "Similarity (higher is better)", y = "Reference Mechanism"
    ) +
    ggplot2::coord_cartesian(xlim = c(min(x$similarity_score, 0), max(x$similarity_score) * 1.1)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5))
}
