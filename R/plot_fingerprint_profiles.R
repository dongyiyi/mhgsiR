#' Plot fingerprint profiles as faceted bars (robust scaling, NA-safe)
#'
#' Visualize one or more fingerprints as faceted bar charts, with a common
#' rescaled axis in [-1, 1]. Missing metrics (NA) are drawn as grey "ghost"
#' bars at 0 and do not affect the scaling. The direction (left/right) always
#' follows the sign of the original value.
#'
#' @param fingerprint_df A data frame returned by `extract_fingerprint()` or
#'   a stacked set of fingerprints (e.g., `reference_fingerprints`).
#'   It should contain the five metric columns (any of the common aliases
#'   are accepted; see below).
#' @param scale_A Optional numeric. If provided, all panels are scaled by the
#'   same amplitude `A` (recommended: compute once from your reference set and
#'   pass it in so case/reference are perfectly comparable). If `NULL`,
#'   the function computes `A` from the supplied data only.
#' @param title Character. Plot title.
#' @param metric_cols Optional character vector of the five preferred metric
#'   column names to draw (defaults are shown below). If `NULL`, the function
#'   will use these defaults and auto-rename common legacy aliases in
#'   `fingerprint_df` to match.
#'
#' Preferred metric names (defaults):
#' \itemize{
#'   \item Abundance–Fitness Slope
#'   \item Stress Response
#'   \item Interaction Strength (R-squared)
#'   \item Diversity–Fitness Link
#'   \item Functional Redundancy
#' }
#' Accepted legacy aliases include e.g. \code{"Abundance_Fitness Slope"},
#' \code{"Interaction Strength (R_squared)"}, \code{"Functional_Redundancy"}, etc.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#'
#' @examples
#' # A: use reference amplitude so panels are comparable
#' A_ref <- max(abs(unlist(reference_fingerprints[
#'   , c("Abundance–Fitness Slope","Stress Response",
#'       "Interaction Strength (R-squared)","Diversity–Fitness Link",
#'       "Functional Redundancy")]
#' )), na.rm = TRUE)
#' plot_fingerprint_profiles(reference_fingerprints, scale_A = A_ref)
#'
#' # B: plot a single case using the same A
#' # plot_fingerprint_profiles(fp_case, scale_A = A_ref)
plot_fingerprint_profiles <- function(
    fingerprint_df,
    scale_A = NULL,
    title   = "Each metric is rescaled to [-1,1] with 0 anchored; NA shown in grey",
    metric_cols = NULL
){
  # ---- 0) Preferred names & alias compatibility (inline, no extra files) ----
  preferred <- c(
    "Abundance_Fitness Slope",
    "Stress Response",
    "Interaction Strength (R_squared)",
    "Diversity_Fitness Link",
    "Functional Redundancy"
  )
  if (is.null(metric_cols)) metric_cols <- preferred

  # map legacy aliases -> preferred
  aliases <- list(
    "Abundance_Fitness Slope"          = c("Abundance_Fitness Slope", "Abundance-Fitness Slope",
                                           "Abundance–Fitness slope"),
    "Stress Response"                  = c("Stress_Response", "Stress response"),
    "Interaction Strength (R_squared)" = c("Interaction Strength (R-squared)",
                                           "Interaction_Strength_R2", "Interaction Strength R2"),
    "Diversity_Fitness Link"           = c("Diversity_Fitness Link", "Diversity–Fitness link"),
    "Functional Redundancy"            = c("Functional_Redundancy")
  )

  df <- fingerprint_df
  nm <- names(df)
  for (pref in preferred) {
    if (!(pref %in% nm)) {
      alt <- aliases[[pref]]
      if (!is.null(alt)) {
        hit <- alt[alt %in% nm]
        if (length(hit) >= 1) {
          names(df)[match(hit[1], names(df))] <- pref
        }
      }
    }
  }
  # any still-missing metric gets an NA column (for ghost bars)
  for (pref in preferred) {
    if (!(pref %in% names(df))) df[[pref]] <- NA_real_
  }

  # facet variable: mechanism or fingerprint id
  facet_col <- if ("mechanism_observed" %in% names(df)) {
    "mechanism_observed"
  } else if ("fingerprint_id" %in% names(df)) {
    "fingerprint_id"
  } else {
    df$fingerprint_id <- paste0("FP_", seq_len(nrow(df)))
    "fingerprint_id"
  }

  # ---- 1) Long format & robust Inf handling ----
  long <- tidyr::pivot_longer(
    data = df,
    cols  = dplyr::all_of(metric_cols),
    names_to = "metric", values_to = "value"
  )

  # replace ±Inf by ±1.1 * finite_max so geoms never blow up
  finite_max <- suppressWarnings(max(long$value[is.finite(long$value)], na.rm = TRUE))
  if (!is.finite(finite_max)) finite_max <- 1
  long <- dplyr::mutate(
    long,
    value = dplyr::case_when(
      is.infinite(value) & value > 0 ~  finite_max * 1.1,
      is.infinite(value) & value < 0 ~ -finite_max * 1.1,
      TRUE ~ value
    )
  )

  # ---- 2) Inline "calc_scale_A": compute amplitude A unless provided ----
  A <- if (is.null(scale_A)) {
    suppressWarnings(max(abs(long$value), na.rm = TRUE))
  } else {
    as.numeric(scale_A)
  }
  if (!is.finite(A) || A == 0) A <- 1

  # ---- 3) Rescale, prepare labels & NA ghost flags ----
  long <- long |>
    dplyr::mutate(
      value_scaled = value / A,
      metric = factor(metric, levels = rev(metric_cols)),
      is_na  = is.na(value),
      label  = ifelse(is_na, "NA", formatC(value, digits = 2, format = "f"))
    )

  # ---- 4) Plot ----
  ggplot2::ggplot(long, ggplot2::aes(x = value_scaled, y = metric)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    # solid bars for non-NA
    ggplot2::geom_col(
      data = subset(long, !is_na),
      ggplot2::aes(fill = value > 0),
      width = 0.7, alpha = 0.9, na.rm = TRUE
    ) +
    # ghost bars for NA at x = 0 (do not affect scaling)
    ggplot2::geom_col(
      data = subset(long, is_na),
      mapping = ggplot2::aes(x = 0, y = metric),
      fill = "grey85", color = "grey70", width = 0.7, alpha = 0.85,
      na.rm = TRUE
    ) +
    # text labels: "NA" for NA; numeric otherwise; placed inside/outside by sign
    ggplot2::geom_text(
      ggplot2::aes(
        label = label,
        x = dplyr::if_else(is.na(value_scaled), 0,
                           dplyr::if_else(value_scaled >= 0,
                                          pmin(value_scaled + 0.04, 1.05),
                                          pmax(value_scaled - 0.04, -1.05))),
        hjust = dplyr::if_else(is.na(value_scaled) | value_scaled >= 0, 0, 1)
      ),
      size = 3.6
    ) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), ncol = 2) +
    ggplot2::scale_fill_manual(values = c(`TRUE`="#0072B2", `FALSE`="#D55E00"), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1.1, 1.1),
                                name = "Rescaled score [-1,1] (original value labeled)") +
    ggplot2::labs(title = title, y = "Diagnostic Metric") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank()
    )
}
