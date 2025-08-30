#' Plot Fingerprint Profiles as Bar Charts (Definitive, Final Version)
#'
#' @details
#' This function visualizes fingerprints as faceted bar charts. It uses a robust
#' global scaling method and now correctly handles `Inf` values, representing them
#' as the maximum on the scaled axis.
#'
#' @param fingerprint_df A data frame from `extract_fingerprint()` or `extract_fingerprint_batch()`.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' plot_fingerprint_profiles(reference_fingerprints)

plot_fingerprint_profiles <- function(fingerprint_df) {

  # --- Input Validation & Preparation ---
  # (This part is correct and remains the same)
  metric_cols <- c(
    "Abundance_Fitness Slope", "Stress Response", "Interaction Strength (R-squared)",
    "Diversity_Fitness Link", "Functional Redundancy"
  )
  if (!all(metric_cols %in% names(fingerprint_df))) {
    stop("Input must contain the 5 core fingerprint metric columns.")
  }
  if ("mechanism_observed" %in% names(fingerprint_df)) {
    group_col <- "mechanism_observed"
  } else if ("fingerprint_id" %in% names(fingerprint_df)) {
    group_col <- "fingerprint_id"
  } else {
    fingerprint_df <- dplyr::mutate(fingerprint_df, fingerprint_id = paste0("FP_", dplyr::row_number()))
    group_col <- "fingerprint_id"
  }

  # --- Data Preparation and Global Scaling ---
  plot_data_long <- fingerprint_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metric_cols),
      names_to = "metric",
      values_to = "value"
    ) %>%
    tidyr::replace_na(list(value = 0))

  # **NEW STEP**: Robustly handle Inf values
  # Find the maximum of all FINITE values in the data
  finite_max <- max(plot_data_long$value[is.finite(plot_data_long$value)], na.rm = TRUE)
  # Replace Inf with a value slightly larger than the finite max, ensuring it will be the top of the scale
  plot_data_long <- plot_data_long %>%
    dplyr::mutate(value = ifelse(is.infinite(value), finite_max * 1.1, value))

  # Step 1: Find the global range across ALL values in the dataset
  global_range <- range(plot_data_long$value, na.rm = TRUE)

  # Step 2: Use this single, global range to scale the entire 'value' column
  plot_data_scaled <- plot_data_long %>%
    dplyr::mutate(
      value_scaled = scales::rescale(value, to = c(-1, 1), from = global_range)
    ) %>%
    dplyr::mutate(metric = factor(metric, levels = rev(metric_cols)))

  # --- Create the Plot ---
  # (The ggplot call remains the same as the previous correct version)
  p <- ggplot2::ggplot(plot_data_scaled, ggplot2::aes(x = value_scaled, y = metric)) +
    ggplot2::geom_col(ggplot2::aes(fill = value > 0), alpha = 0.9, width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(value, 2),
          x = value_scaled + ifelse(value > 0, 0.02, -0.02),
          hjust = ifelse(value > 0, 0, 1)),
      size = 3.5
    ) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", group_col)), ncol = 2) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(-1.1, 1.1), name = "Rescaled Metric Score (Original Value Labeled)") +
    ggplot2::labs(
      title = "Comparison of Symbiotic Mechanism Fingerprints",
      subtitle = "Each metric is rescaled to a common [-1, 1] axis for comparison.",
      y = "Diagnostic Metric"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30")
    )

  return(p)
}
