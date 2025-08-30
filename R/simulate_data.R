#' @importFrom magrittr %>%
NULL

#' Simulate Symbiotic Interaction Data
#'
#' This is the core simulation engine for the mhgsiR package. It generates a
#' dataset of host fitness based on symbiont abundance under one of four
#' core mechanistic archetypes.
#'
#' @details
#' (The excellent philosophy text you wrote remains here...)
#'
#' @param mechanism A character string specifying the symbiotic mechanism.
#' @param n A single positive integer for the number of hosts to simulate.
#' @param noise A single non-negative numeric value for the noise standard deviation.
#' @param is_stress A logical value (TRUE/FALSE) for the stress condition.
#' @param verbose A logical value. If TRUE, prints simulation parameters.
#' @param seed An optional integer to set the random seed for reproducibility.
#' @param custom_params An optional list with `intercept` and `slope`.
#' @param add_label A logical value. If TRUE, adds a human-readable, ordered factor column.
#' @param return_invisible A logical value. If TRUE, the function returns the
#'   data frame invisibly, preventing auto-printing in the console.
#'
#' @return A tibble with host-symbiont data. The ground truth parameters of the
#'   simulation are stored as an attribute named "params".
#' @export
#' @examples
#' # Simulate and pipe directly to a plot without printing the data
#' simulate_data(mechanism = "Defensive", seed = 42, return_invisible = TRUE) %>%
#'   ggplot2::ggplot(aes(x = abundance, y = fitness)) +
#'   ggplot2::geom_point()

simulate_data <- function(mechanism = c("Nutritional", "Defensive", "Neutral", "FreeRider", "Custom"),
                          n = 150,
                          noise = 0.15,
                          is_stress = FALSE,
                          verbose = TRUE,
                          seed = NULL,
                          custom_params = NULL,
                          add_label = FALSE,
                          include_id = TRUE,
                          return_invisible = FALSE) {

  # --- Input Validation ---
  mechanism <- match.arg(mechanism)
  # ... (All previous input validation checks remain here) ...
  if (mechanism == "Custom" && is.null(custom_params)) {
    stop("You must supply `custom_params` when using mechanism = 'Custom'.")
  }
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n %% 1 != 0) {
    stop("`n` must be a single positive integer.")
  }
  if (!is.numeric(noise) || length(noise) != 1 || noise < 0) {
    stop("`noise` must be a single non-negative numeric value.")
  }


  # --- Seed for Reproducibility ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Parameter Setup ---
  if (!is.null(custom_params)) {
    # ... (logic for custom_params) ...
    params <- custom_params
    mechanism <- "Custom"
  } else {
    # ... (default parameter logic) ...
    params <- dplyr::case_when(
      mechanism == "Nutritional" & is_stress  ~ list(intercept = 0.3, slope = 0.6),
      mechanism == "Nutritional" & !is_stress ~ list(intercept = 0.5, slope = 0.0),
      mechanism == "Defensive" & is_stress  ~ list(intercept = 0.2, slope = 0.7),
      mechanism == "Defensive" & !is_stress ~ list(intercept = 0.5, slope = -0.1),
      mechanism == "Neutral"                 ~ list(intercept = 0.5, slope = 0.0),
      mechanism == "FreeRider"                ~ list(intercept = 0.5, slope = -0.3)
    )
  }

  # ... (verbose logging) ...
  if (verbose) {
    condition_str <- ifelse(is_stress, "Stress", "No Stress")

    # First, construct the entire message string
    story <- sprintf(
      "Simulating: '%s' mechanism | Condition: '%s'\n  - Parameters: intercept=%.2f, slope=%.2f",
      mechanism,
      condition_str,
      params$intercept,
      params$slope
    )

    # Then, print the clean string with message()
    message(story)
  }

  # --- Data Generation ---
  df <- tibble::tibble(
    abundance = stats::runif(n, 0, 1)
  )

  if (include_id) {
    df <- dplyr::mutate(df, host_id = 1:n) %>%
      dplyr::relocate(host_id, .before = abundance)
  }

  df <- df %>%
    dplyr::mutate(
      fitness = params$intercept + params$slope * abundance + stats::rnorm(n, 0, noise),
      mechanism = mechanism,
      condition = ifelse(is_stress, "Stress", "No Stress")
    )

  # --- Add Optional Ordered Factor Label ---
  if (add_label) {
    label_text <- switch(mechanism,
                         Nutritional = "Nutritional (Benefit under stress)",
                         Defensive = "Defensive (Cost/Benefit)",
                         Neutral = "Neutral (No effect)",
                         FreeRider = "Free-Rider (Cost only)",
                         Custom = "Custom (User-defined)"
    )

    # Define the desired order for the factor levels
    label_levels <- c(
      "Nutritional (Benefit under stress)",
      "Defensive (Cost/Benefit)",
      "Neutral (No effect)",
      "Free-Rider (Cost only)",
      "Custom (User-defined)"
    )

    df <- dplyr::mutate(df, mechanism_label = factor(label_text, levels = label_levels))
  }

  # --- Embed Parameters as Metadata Attribute ---
  attr(df, "params") <- list(
    mechanism = mechanism,
    intercept = params$intercept,
    slope = params$slope,
    stress = is_stress,
    noise = noise,
    n = n,
    seed = seed
  )

  # --- Return Output ---
  if (return_invisible) {
    return(invisible(df))
  } else {
    return(df)
  }
}
