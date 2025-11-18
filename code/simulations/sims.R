library(tidyverse)
library(here)
library(future)
library(furrr)
library(progressr)

source(here("code/simulations-new/sim_helpers.R"))

# Full set of simulations ------------------------------------------------------
# %%
set.seed(20241003)
sims <- expand_grid(
  expand_grid(
    N = c(50, 300),
    T = c(5, 10, 15)
  ),
  expand_grid(
    tibble(
      dgp_number = c(1, 2, 3),
      avg_diff_in_loading = c(0, -0.5, -0.5),
      eta_y = c(0, 0, 0),
      tau_x = c(0, 0, 1),
    ),
    eps_autocorrelated = c(FALSE, TRUE),
    fact_label = c("Constant", "Linear Trend", "Autocorrelated")
  )
)

# sims <- expand_grid(
#   expand_grid(
#     N = c(50, 100, 300),
#     T = c(5, 10, 15)
#   ),
#   tibble(
#     dgp_num = 1:3,
#     fact_label = rep(c("Constant", "Linear Trend", "Autocorrelated"), each = 1),
#     avg_diff_in_loading = rep(c(0), times = 3),
#     eta_y = rep(c(0), times = 3),
#     tau_x = rep(c(0), times = 3)
#   )
# )

# %%
set.seed(20241114)
future::plan("multicore", workers = 6)
nsim <- 250

sim_results <- seq_len(nrow(sims)) |>
  # purrr::map(
  furrr::future_map(
    function(i) {
      sim <- sims[i, ]
      # catma("sim: {i}\n")
      res <- run_simulations(
        nsim = nsim,
        N = sim$N,
        T = sim$T,
        dgp_number = sim$dgp_number,
        fact_label = sim$fact_label,
        avg_diff_in_loading = sim$avg_diff_in_loading,
        eta_y = sim$eta_y,
        tau_x = sim$tau_x,
        eps_autocorrelated = sim$eps_autocorrelated
      )

      return(res)
    },
    .options = furrr::furrr_options(seed = TRUE)
  ) |>
  list_rbind() |>
  as_tibble()

sim_results <- sim_results |>
  mutate(
    fact_label = fact_label |>
      as.character() |>
      fct(levels = c("Constant", "Linear Trend", "Autocorrelated")),
    err_label = ifelse(eps_autocorrelated, "AR(1)", "Independent"),
    dgp_label = sprintf(
      "$\\kappa = %s$, $\\tau_x = %s$",
      avg_diff_in_loading,
      tau_x
    )
  )


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs::dir_create(here("data/simulations-new/"))
write_csv(
  x = sim_results,
  here("data/simulations-new/sim_results.csv")
)
