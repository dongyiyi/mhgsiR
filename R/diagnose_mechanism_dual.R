# This internal helper contains the core calculation logic.
# It is NOT exported, so users will not see it.
.diagnose_core <- function(observed_fingerprint,
                           reference_library = mhgsiR::reference_fingerprints,
                           method = c("cosine", "correlation", "euclidean"),
                           correlation_type = c("pearson","spearman","kendall"),
                           allow_rescale = FALSE,
                           na_strategy = c("drop_metric", "error")) {

  method <- match.arg(method)
  correlation_type <- match.arg(correlation_type)
  na_strategy <- match.arg(na_strategy)
  metric_names <- .mhgsi_metric_names

  # --- (B) Mechanism 兼容处理 ---
  if (!"Mechanism" %in% names(reference_library)) {
    alt <- c("mechanism_observed","mechanism","mechanism_name")
    hit <- alt[alt %in% names(reference_library)][1]
    if (!is.na(hit)) {
      reference_library <- dplyr::rename(reference_library, Mechanism = !!rlang::sym(hit))
    } else {
      stop("`reference_library` must contain a `Mechanism` (or legacy) column.", call. = FALSE)
    }
  }
  # 输入校验
  if (!all(metric_names %in% names(observed_fingerprint)) ||
      !all(metric_names %in% names(reference_library))) {
    stop("Inputs must contain all 5 core metric columns.", call. = FALSE)
  }

  # 数据准备
  observed_vec     <- as.numeric(observed_fingerprint[1, metric_names])
  reference_matrix <- as.matrix(reference_library[, metric_names, drop = FALSE])

  # --- 更稳健：把 NA/NaN/Inf 都视作“不可用” ---
  is_finite_obs     <- is.finite(observed_vec)
  is_finite_ref_col <- apply(is.finite(reference_matrix), 2, all)

  if (na_strategy == "drop_metric") {
    keep_cols <- is_finite_obs & is_finite_ref_col
    if (!any(keep_cols)) stop("All metrics are NA/NaN/Inf; cannot compute.", call. = FALSE)
    observed_vec_clean     <- observed_vec[keep_cols]
    reference_matrix_clean <- reference_matrix[, keep_cols, drop = FALSE]
    used_metrics           <- metric_names[keep_cols]
  } else {
    if (!all(is_finite_obs) || !all(is_finite_ref_col))
      stop("Non-finite values detected. Use na_strategy='drop_metric'.", call. = FALSE)
    observed_vec_clean     <- observed_vec
    reference_matrix_clean <- reference_matrix
    used_metrics           <- metric_names
  }

  # 可选 z-score
  if (isTRUE(allow_rescale)) {
    cm <- rbind(observed_vec_clean, reference_matrix_clean)
    sm <- scale(cm); sm[is.nan(sm)] <- 0
    observed_vec_clean     <- sm[1, ]
    reference_matrix_clean <- sm[-1, , drop = FALSE]
  }

  # 相似度
  cosine_sim <- function(a,b){
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (d < .Machine$double.eps) return(0)
    sum(a*b)/d
  }
  similarity_fn <- switch(method,
                          correlation = function(a,b){
                            if (stats::sd(a)==0 || stats::sd(b)==0) return(0)
                            suppressWarnings(stats::cor(a, b, method = correlation_type, use = "pairwise.complete.obs"))
                          },
                          euclidean   = function(a,b) -sqrt(sum((a-b)^2)),
                          cosine      = cosine_sim
  )

  scores <- apply(reference_matrix_clean, 1, function(r) similarity_fn(observed_vec_clean, r))

  similarity_metric <- if (identical(method, "correlation")) paste0("cor_", correlation_type) else method

  result <- tibble::tibble(
    Mechanism         = reference_library$Mechanism,
    similarity_score  = as.numeric(scores),
    method            = method,
    correlation_type = if (identical(method,"correlation")) correlation_type else NA_character_,
    similarity_metric= if (identical(method,"correlation")) paste0("cor_", correlation_type) else method,
    n_metrics         = length(used_metrics)
  ) |>
    dplyr::arrange(dplyr::desc(similarity_score)) |>
    dplyr::mutate(
      rank     = dplyr::row_number(),
      tie_rank = dplyr::min_rank(dplyr::desc(similarity_score))
    )

  attr(result, "used_metrics")    <- used_metrics
  attr(result, "observed_vector") <- as.numeric(observed_vec_clean)
  attr(result, "n_metrics")       <- length(used_metrics)
  attr(result, "similarity_metric") <- similarity_metric
  result
}


#' Diagnose a Symbiotic Mechanism with Automated Method Selection
#'
#' @description
#' The primary diagnostic function of the mhgsiR package. It compares an
#' observed fingerprint to a reference library and provides a ranked diagnosis.
#' The `method = "auto"` setting performs a robustness check across multiple
#' similarity metrics and recommends the most suitable one.
#'
#' @details
#' (The excellent "Choice Criteria" text you wrote will be placed here in the final documentation...)
#'
#' @param observed_fingerprint A 1-row tibble from `extract_fingerprint()`.
#' @param reference_library A data.frame of reference fingerprints. Defaults to `mhgsiR::reference_fingerprints`.
#' @param method The similarity metric. Defaults to `"auto"`, which performs a
#'   concordance analysis and selects the best method. Other options are
#'   `"cosine"`, `"correlation"`, and `"euclidean"`.
#' @param allow_rescale Logical. If FALSE (default), compares raw values. If TRUE,
#'   applies z-score scaling. This setting is used by the chosen method.
#'
#' @return A list of class `mhgsi_diagnosis` containing:
#'   \item{summary}{A ranked tibble with the diagnostic results from the recommended method.}
#'   \item{concordance}{A summary of the auto-selection process, including the
#'     recommended method, reason, and agreement between all methods.}
#'   \item{all_results}{A list containing the full diagnosis tables for each method tested.}
#' @export
#'
#' @examples
#' obs_fp <- extract_fingerprint(simulate_data("Defensive", is_stress = TRUE, seed = 101))
#' diagnosis <- diagnose_mechanism(obs_fp)
#' print(diagnosis)
#' plot(diagnosis)

diagnose_mechanism_dual <- function(observed_fingerprint,
                                    reference_library = mhgsiR::reference_fingerprints,
                                    method = c("auto", "cosine", "correlation", "euclidean"),
                                    correlation_type = c("pearson","spearman","kendall"),
                                    allow_rescale = FALSE,
                                    report = c("quiet","message","verbose")) {

  method <- match.arg(method)
  correlation_type <- match.arg(correlation_type)
  report <- match.arg(report)

  # -------- manual path
  if (!identical(method, "auto")) {
    res <- .diagnose_core(
      observed_fingerprint, reference_library,
      method = method,
      correlation_type = correlation_type,
      allow_rescale = allow_rescale
    )

    out <- list(
      summary = res,
      concordance = list(
        recommended_method = method,
        correlation_type   = if (identical(method, "correlation")) correlation_type else NA_character_,
        reason             = "Method manually specified by user.",
        top1_agreement     = TRUE,
        top1_mechanisms    = res$Mechanism[1],
        rank_correlation_matrix = matrix(1, nrow = 1, ncol = 1,
                                         dimnames = list(attr(res,"similarity_metric"),
                                                         attr(res,"similarity_metric"))),
        n_metrics_used     = attr(res, "n_metrics"),
        used_metrics       = attr(res, "used_metrics")
      ),
      all_results = setNames(list(res), attr(res,"similarity_metric"))
    )
    class(out) <- c("mhgsi_diagnosis_dual", "mhgsi_diagnosis")

    if (report != "quiet") {
      msg <- sprintf("method = '%s'%s • allow_rescale = %s",
                     method,
                     if (identical(method,"correlation")) paste0(" (", correlation_type, ")") else "",
                     allow_rescale)
      message(msg)
    }
    return(out)
  }

  # -------- AUTO path: evaluate multiple metrics (with rescale=TRUE for fair comparison)
  cand <- list(
    cor_pearson = .diagnose_core(observed_fingerprint, reference_library,
                                 method="correlation", correlation_type="pearson", allow_rescale=FALSE),
    cor_spearman = .diagnose_core(observed_fingerprint, reference_library,
                                  method="correlation", correlation_type="spearman", allow_rescale=FALSE),
    cosine = .diagnose_core(observed_fingerprint, reference_library,
                            method="cosine", allow_rescale=FALSE),
    euclidean = .diagnose_core(observed_fingerprint, reference_library,
                               method="euclidean", allow_rescale=FALSE)
  )

  # rank concordance
  rank_tbl <- Reduce(function(a,b) dplyr::full_join(a,b,by="Mechanism"), lapply(names(cand), function(nm) {
    dplyr::select(cand[[nm]], Mechanism, !!rlang::sym(paste0(nm, "_rank")) := rank)
  }))
  rank_mat <- as.matrix(rank_tbl[ , grep("_rank$", names(rank_tbl), value = TRUE)])
  colnames(rank_mat) <- names(cand)
  rank_corr <- suppressWarnings(stats::cor(rank_mat, method="spearman", use="pairwise.complete.obs"))

  # top-1 mechanisms and margins
  top1_mechanisms <- vapply(cand, function(dx) dx$Mechanism[1], character(1))
  top1_agreement  <- length(unique(top1_mechanisms)) == 1
  top1_margin <- vapply(cand, function(dx) {
    if (nrow(dx) < 2) return(NA_real_)
    dx$similarity_score[1] - dx$similarity_score[2]
  }, numeric(1))

  # average agreement of each candidate with others
  avg_agreement <- vapply(seq_along(cand), function(i) {
    mean(rank_corr[i, -i], na.rm = TRUE)
  }, numeric(1))
  nmet <- attr(cand[[1]], "n_metrics")

  .choose_corr <- function(avg_agreement, default = "cor_spearman") {
    a <- unname(avg_agreement["cor_spearman"])
    b <- unname(avg_agreement["cor_pearson"])
    if (is.na(a) && is.na(b)) return(default)
    if (is.na(a)) return("cor_pearson")
    if (is.na(b)) return("cor_spearman")
    if (a > b) "cor_spearman" else if (a < b) "cor_pearson" else default
  }

  # 用 -Inf 代替 NA，保证 max/排序能工作（不覆盖原值时可复制向量）
  avg_agreement_f <- avg_agreement
  avg_agreement_f[is.na(avg_agreement_f)] <- -Inf

  pick_name   <- NULL
  pick_reason <- NULL

  if (top1_agreement) {
    cor_choice <- .choose_corr(avg_agreement)
    pick_name  <- cor_choice
    pick_reason <- "All methods agreed on Top-1; chose correlation with higher rank agreement (NA-safe; tie→spearman)."

  } else if (nmet <= 3) {
    pick_name  <- "cosine"
    pick_reason <- "Methods disagreed with ≤3 metrics; cosine is more robust for low-dimensional fingerprints."

  } else {
    # 先在有限值里挑最高一致性的候选；若全是 -Inf（原本是 NA），退回全体
    candidates <- names(avg_agreement_f)[is.finite(avg_agreement_f)]
    if (length(candidates) == 0) candidates <- names(avg_agreement_f)

    tied <- candidates[avg_agreement_f[candidates] == max(avg_agreement_f[candidates], na.rm = TRUE)]

    if (length(tied) > 1) {
      # 用 Top-1 margin 打破平手；仍平手则优先 spearman
      margins <- top1_margin[tied]
      tied <- tied[margins == max(margins, na.rm = TRUE)]
      if (length(tied) > 1 && "cor_spearman" %in% tied) tied <- "cor_spearman"
    }
    pick_name <- tied[1]
    pick_reason <- "Methods disagreed; picked highest average rank agreement (NA-safe; tie→larger Top-1 margin, then spearman)."
  }

  # 映射为最终 method / correlation_type
  if (pick_name %in% c("cor_pearson","cor_spearman")) {
    pick_method <- "correlation"
    pick_corr   <- if (pick_name == "cor_pearson") "pearson" else "spearman"
  } else if (pick_name == "cosine") {
    pick_method <- "cosine";    pick_corr <- NA_character_
  } else {
    pick_method <- "euclidean"; pick_corr <- NA_character_
  }

  # Final pass with user's allow_rescale
  final <- .diagnose_core(
    observed_fingerprint, reference_library,
    method = pick_method,
    correlation_type = if (is.na(pick_corr)) "pearson" else pick_corr,
    allow_rescale = allow_rescale
  )

  # assemble concordance / all_results
  out <- list(
    summary = final,
    concordance = list(
      recommended_method      = pick_method,
      correlation_type        = if (identical(pick_method,"correlation")) pick_corr else NA_character_,
      reason                  = pick_reason,
      top1_agreement          = top1_agreement,
      top1_mechanisms         = top1_mechanisms,
      method_scores           = tibble::tibble(
        method = names(cand),
        avg_rank_agreement = unname(avg_agreement),
        top1_margin = unname(top1_margin),
        top1 = unname(top1_mechanisms)
      ) |> dplyr::arrange(dplyr::desc(avg_rank_agreement), dplyr::desc(top1_margin)),
      rank_correlation_matrix = rank_corr,
      n_metrics_used          = nmet,
      used_metrics            = attr(cand[[1]], "used_metrics")
    ),
    all_results = cand
  )
  class(out) <- c("mhgsi_diagnosis_dual", "mhgsi_diagnosis")

  if (report != "quiet") {
    pretty_metric <- if (identical(pick_method,"correlation")) paste0("cor_", pick_corr) else pick_method
    msg <- sprintf("method = 'auto' → using %s • allow_rescale = %s", pretty_metric, allow_rescale)
    message(msg)
    if (report == "verbose") {
      utils::capture.output(print(out$concordance$method_scores))
    }
  }

  out
}


#' @export
#' @method plot mhgsi_diagnosis_dual
plot.mhgsi_diagnosis_dual <- function(x, ...) {
  res <- x$summary
  sublab <- sprintf(
    "Recommended: %s%s  •  %s",
    x$concordance$recommended_method,
    if (!is.na(x$concordance$correlation_type)) paste0(" (", x$concordance$correlation_type, ")") else "",
    x$concordance$reason
  )
  ggplot2::ggplot(res, ggplot2::aes(x = similarity_score, y = stats::reorder(Mechanism, similarity_score))) +
    ggplot2::geom_col(ggplot2::aes(fill = rank == 1L), alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", similarity_score)), hjust = -0.15, size = 3.5) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey70"), guide = "none") +
    ggplot2::labs(
      title = "Mechanism Diagnosis",
      subtitle = sublab,
      x = "Similarity (higher is better)", y = "Mechanism"
    ) +
    ggplot2::theme_classic(base_size = 12)
}
