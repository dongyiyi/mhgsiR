#' Extract Fingerprints for Theoretical Mechanisms
#'
#' Generates a standardized library of mechanistic fingerprints. By default, it
#' returns the idealized, noise-free theoretical fingerprints.
#'
#' @param mechanisms A character vector of mechanisms to generate.
#' @param idealized Logical. If TRUE (the default), returns the perfect,
#'   noise-free theoretical fingerprints. If FALSE, runs a noisy simulation
#'   to generate the fingerprints.
#' @param n Number of hosts (only used if `idealized = FALSE`).
#' @param noise Noise level (only used if `idealized = FALSE`).
#' @param seed Optional seed for reproducibility (only used if `idealized = FALSE`).
#'
#' @return A tibble of theoretical fingerprints.
#' @export

extract_fingerprint_batch <- function(mechanisms = c("Nutritional", "Defensive", "Neutral", "FreeRider"),
                                                idealized = TRUE,
                                                n = 10000, # Use large n for stable simulations
                                                noise = 0.1, # Use low noise
                                                seed = 123) {

  if (idealized) {
    # --- IDEALIZED MODE: Return the perfect theoretical values ---

    # These are the direct mathematical consequences of our model's parameters
    reference_library <- tibble::tribble(
      ~fingerprint_id, ~`Abundance_Fitness Slope`, ~`Stress Response`, ~`Interaction Strength (R-squared)`, ~`Diversity_Fitness Link`, ~`Functional Redundancy`, ~mechanism_observed,
      "ref_Nutritional",                      0.0,               0.6,                                 1.0,                     NA,                      NA,       "Nutritional",
      "ref_Defensive",                       -0.1,               0.8,                                 1.0,                     NA,                      NA,         "Defensive",
      "ref_Neutral",                          0.0,               0.0,                                 0.0,                     NA,                      NA,           "Neutral",
      "ref_FreeRider",                       -0.3,               0.0,                                 1.0,                     NA,                      NA,         "FreeRider"
    )
    reference_library <- reference_library %>%
      dplyr::filter(mechanism_observed %in% mechanisms)


  } else {
    # --- SIMULATION MODE with ROBUST SEEDING ---

    # We use purrr::imap_dfr which provides both the value (.x) and the index (.y)
    reference_library <- purrr::imap_dfr(mechanisms, function(mech, index) {

      # Generate a unique, deterministic seed for EACH mechanism
      mechanism_seed <- if (!is.null(seed)) seed + index else NULL

      message("Running high-n simulation for: ", mech, " (using seed ", mechanism_seed, ")")

      # Simulate data for both conditions using the specific seed for this mechanism
      combined_data <- dplyr::bind_rows(
        simulate_data(mechanism = mech, is_stress = TRUE, n = n, noise = noise, seed = mechanism_seed, verbose = FALSE),
        simulate_data(mechanism = mech, is_stress = FALSE, n = n, noise = noise, seed = mechanism_seed, verbose = FALSE)
      )

      # Extract the single fingerprint from this combined dataset
      extract_fingerprint(combined_data, fingerprint_id = paste0("sim_", mech))
    })
  }

  return(reference_library)
}
