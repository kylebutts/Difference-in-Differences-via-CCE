#' # Main script to run the complete project
#'
#' `render_file` will run the R file and log the results of the script
#' in the `logbook`. The logbook can be viewed with
#' `quarto preview logbook` from the terminal.
library(here)
source(here("logbook/render_file.R"))

# Trade Liberalization ---------------------------------------------------------
1 && render_file(here(
  "code/Trade-Liberalization-and-Markup-Dispersion/analysis.R"
))
1 && render_file(here(
  "code/Trade-Liberalization-and-Markup-Dispersion/analysis_interpolate_2003.R"
))

# Simulations ------------------------------------------------------------------
0 && render_file(here("code/simulations/misspecified_exposure_mapping.R"))
