# mhgsiR
*From Stories to Signatures: A Framework for Inferring Symbiotic Processes*

<!-- badges: start -->
[![R-CMD-check](https://github.com/dongyiyi/mhgsiR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dongyiyi/mhgsiR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`mhgsiR` is an R package for **simulating, extracting, and diagnosing macro-fingerprints** of host–symbiont interactions.  
It translates ecological *stories* into standardized, diagnosable **signatures**, enabling mechanism inference, comparison, and robustness analysis.

## Features
- **Simulate** host–symbiont datasets under four mechanistic archetypes: Nutritional, Defensive, FreeRider, Neutral  
- **Extract** standardized 5-dimensional macro-fingerprints from observational or simulated data  
- **Diagnose** mechanisms by comparing observed fingerprints to a theoretical reference library  
- **Visualize** conceptual scenarios, fingerprint profiles, and diagnostic bar plots  
- **Assess robustness** via power analysis (sample size × noise maps)

## Installation
```r
# install.packages("devtools")
devtools::install_github("dongyiyi/mhgsiR")
