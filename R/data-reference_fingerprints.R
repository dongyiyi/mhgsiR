#' Reference Fingerprints for Mechanistic Archetypes
#'
#' Idealized, noise-free macro-fingerprint profiles for four archetypes
#' (Defensive, Nutritional, FreeRider, Neutral). Used as the default
#' comparison library in \code{diagnose_mechanism_dual()} and related functions.
#'
#' @format A data frame with 4 rows and 7 columns:
#' \describe{
#'   \item{fingerprint_id}{Character ID (e.g., "ref_Defensive").}
#'   \item{Mechanism}{Mechanism name.}
#'   \item{abundance_fitness_slope}{Numeric.}
#'   \item{stress_response}{Numeric.}
#'   \item{interaction_strength_r2}{Numeric.}
#'   \item{diversity_fitness_link}{Numeric.}
#'   \item{functional_redundancy}{Numeric.}
#' }
#' @source Constructed internally for the mhgsiR package.
"reference_fingerprints"
