#' Extract a Multi-dimensional Mechanistic Fingerprint
#'
#' From a host-symbiont dataset, this function calculates the five core
#' diagnostic metrics of the MHG-SI framework.
#'
#' @details
#' This is the analytical core of the mhgsiR package. It synthesizes raw data
#' into the standardized 'macro-fingerprint' vector required for diagnosis.
#' Note that `Diversity–Fitness Link` and `Functional Redundancy` are
#' currently returned as `NA` placeholders.
#'
#' @param df A tibble with columns: `abundance`, `fitness`, and `condition`.
#' @param return_full Logical. If TRUE, returns a detailed tibble with all
#'   intermediate metrics calculated for each condition. Defaults to FALSE.
#' @param fingerprint_id An optional character string to label the output. If NULL,
#'   a unique ID will be generated based on the system time.
#'
#' @return A single-row tibble containing the complete 5D fingerprint (default),
#'   or a multi-row tibble with detailed metrics if `return_full = TRUE`.
#' @export
#'
#' @examples
#' combined_data <- dplyr::bind_rows(
#'   simulate_data(mechanism = "Defensive", is_stress = TRUE, seed = 42),
#'   simulate_data(mechanism = "Defensive", is_stress = FALSE, seed = 42)
#' )
#' extract_fingerprint(combined_data)

extract_fingerprint <- function(df, fingerprint_id = NULL) {
  # --- Input Checks ---
  required_cols <- c("abundance", "fitness")
  if (!all(required_cols %in% names(df))) {
    stop("Input `df` must contain 'abundance' and 'fitness' columns.")
  }
  if (!"condition" %in% names(df)) {
    df$condition <- "No Stress"
  }

  # Helper function remains the same
  extract_lm_stats <- function(sub_df) {
    if (nrow(sub_df) < 3) return(tibble::tibble(slope = NA, r_squared = NA))
    fit <- stats::lm(fitness ~ abundance, data = sub_df)
    summary_fit <- summary(fit)
    tibble::tibble(slope = unname(stats::coef(fit)[2]), r_squared = summary_fit$r.squared)
  }

  full_metrics <- df %>%
    dplyr::group_by(condition) %>%
    tidyr::nest() %>%
    dplyr::mutate(stats = purrr::map(data, extract_lm_stats)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(stats)

  # --- Synthesize the 5D Fingerprint ---
  slope_stress <- (full_metrics %>% dplyr::filter(condition == "Stress"))$slope
  if (length(slope_stress) == 0) slope_stress <- 0

  slope_no_stress <- (full_metrics %>% dplyr::filter(condition == "No Stress"))$slope
  if (length(slope_no_stress) == 0) slope_no_stress <- 0

  # **CORRECTION 1: Extract R-squared, not variance**
  # Extract R-squared from BOTH conditions
  r_sq_stress <- (full_metrics %>% dplyr::filter(condition == "Stress"))$r_squared
  if (length(r_sq_stress) == 0) r_sq_stress <- 0
  r_sq_no_stress <- (full_metrics %>% dplyr::filter(condition == "No Stress"))$r_squared
  if (length(r_sq_no_stress) == 0) r_sq_no_stress <- 0

  # Define the metric as the MAXIMUM of the two, reflecting the strongest interaction
  metric_interaction_strength <- max(r_sq_stress, r_sq_no_stress)


  # (Robustly check for a single, unique mechanism name)
  mechanism_observed <- if ("mechanism" %in% names(df)) {
    m <- unique(df$mechanism)
    if (length(m) == 1) m else NA_character_
  } else {
    NA_character_
  }

  # (Generate fingerprint_id if not provided)
  if (is.null(fingerprint_id)) {
    fingerprint_id <- paste0("fp_", as.integer(Sys.time()))
  }

  # Assemble the final, single-row fingerprint
  fingerprint <- tibble::tibble(
    fingerprint_id = fingerprint_id,
    `Abundance_Fitness Slope` = slope_no_stress,
    `Stress Response` = slope_stress - slope_no_stress,
    `Interaction Strength (R_squared)` = metric_interaction_strength, # Use the new, correct value
    `Diversity_Fitness Link` = NA_real_,
    `Functional Redundancy` = NA_real_,
    mechanism_observed = mechanism_observed
  ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

  return(fingerprint)
}
