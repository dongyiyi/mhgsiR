#' Plot Simulated Symbiotic Scenarios
#'
#' Visualizes simulated host-symbiont data under various mechanistic archetypes
#' and stress conditions, producing a clean, publication-ready multi-panel plot.
#'
#' @details
#' This function serves as a powerful visualization tool to compare the conceptual
#' fingerprints of different symbiotic mechanisms. The `facet_by = 'wrap'` option
#' provides a compact view, while `facet_by = 'grid'` creates a stricter
#' matrix comparing every mechanism across both stress conditions.
#'
#' @param mechanisms A character vector of symbiotic mechanisms to simulate. Defaults to all four.
#' @param n Number of hosts to simulate per mechanism × condition (default = 100).
#' @param noise Standard deviation of noise to add to host fitness (default = 0.15).
#' @param seed Optional integer for reproducibility.
#' @param facet_by One of "wrap" or "grid". "wrap" facets by mechanism only. "grid" facets by condition × mechanism.
#' @param return_data Logical. If TRUE, returns the full simulated dataset instead of the plot. Default = FALSE.
#'
#' @return A ggplot object (or a tibble if return_data = TRUE).
#' @export
#'
#' @examples
#' # Basic usage to see all scenarios
#' plot_scenarios()
#'
#' # Use a grid layout for a more structured comparison
#' plot_scenarios(facet_by = "grid")
#'
#' # Generate the data for custom analysis instead of plotting
#' scenario_df <- plot_scenarios(return_data = TRUE)

plot_scenarios <- function(mechanisms = c("Nutritional", "Defensive", "Neutral", "FreeRider"),
                           n = 100,
                           noise = 0.15,
                           seed = NULL,
                           facet_by = c("wrap", "grid"),
                           return_data = FALSE) {

  facet_by <- match.arg(facet_by)

  # --- Set seed once for reproducibility ---
  if (!is.null(seed)) set.seed(seed)

  # --- Define condition-color palette for consistency ---
  mhgsi_colors <- c("No Stress" = "#999999", "Stress" = "#D55E00")

  # --- Create simulation grid and simulate all scenarios ---
  scenario_grid <- tidyr::expand_grid(
    mechanism = mechanisms,
    is_stress = c(FALSE, TRUE)
  )

  # A safe version of our simulation function
  safe_simulate <- purrr::possibly(simulate_data, otherwise = NULL)

  sim_data <- purrr::pmap_dfr(
    scenario_grid,
    safe_simulate,
    n = n,
    noise = noise,
    add_label = TRUE,
    verbose = FALSE,
    return_invisible = TRUE
  ) |>
    # Filter out any potential simulation failures
    dplyr::filter(!is.null(mechanism))

  # --- Optional return of raw data ---
  if (return_data) return(sim_data)

  # --- Create base plot ---
  p <- ggplot2::ggplot(sim_data, ggplot2::aes(x = abundance, y = fitness, color = condition)) +
    ggplot2::geom_point(alpha = 0.4, size = 1) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
    ggplot2::scale_color_manual(values = mhgsi_colors, name = "Condition:") +
    ggplot2::labs(
      title = "Comparison of Symbiotic Mechanism Scenarios",
      subtitle = "Each mechanism simulated under Stress and No Stress conditions",
      x = "Symbiont Abundance",
      y = "Host Fitness"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      legend.position = "bottom"
    )

  # --- Apply chosen facet layout ---
  if (facet_by == "wrap") {
    p <- p + ggplot2::facet_wrap(~ mechanism_label, nrow = 2)
  } else {
    p <- p + ggplot2::facet_grid(condition ~ mechanism_label)
  }

  return(p)
}
