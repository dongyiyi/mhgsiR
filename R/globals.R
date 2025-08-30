# In R/globals.R

#' Internal vector of official metric names
#'
#' This vector defines the official column names for the five macro-fingerprint
#' metrics used throughout the mhgsiR package. Centralizing the names here
#' ensures consistency across functions and makes future renaming trivial.
#'
#' @keywords internal
.mhgsi_metric_names <- c(
  "Abundance–Fitness Slope",
  "Stress Response",
  "Interaction Strength (R-squared)",
  "Diversity–Fitness Link",
  "Functional Redundancy"
)

utils::globalVariables(c(
  "similarity_score","Mechanism","mechanism","mechanism_observed",
  "avg_rank_agreement","condition","data","stats","where",
  "value","metric","value_scaled","abundance","fitness",
  "host_id"
))
