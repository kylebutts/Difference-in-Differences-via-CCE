# %%
library(tidyverse)
library(tinytable)
library(here)
library(kfbmisc) # remotes::install_github("kylebutts/kfbmisc")

# fs::dir_create(here("out/tables/simulation/"))

# sim_results <- here("data/simulations/sim_results.csv") |>
#   read_csv(show_col_types = FALSE)

sim_results <- here("data/simulations-new/sim_results.csv") |>
  read_csv(show_col_types = FALSE) |>
  filter(estimator != "gsynth + unit FE") |> # Don't use for now
  mutate(
    fact_label = fact_label |>
      fct(levels = c("Constant", "Linear Trend", "Autocorrelated")),
    dgp_str = sprintf(
      "$\\kappa = %s$, $\\tau_x = %s$",
      avg_diff_in_loading,
      tau_x
    ),
  )

# res <- sim_results |>
#   filter(
#     fact_label == "Linear Trend",
#     avg_diff_in_loading == -0.5,
#     tau_x == 1
#   )

make_table <- function(res, tbl_caption) {
  table_order <- c(
    "TECCE" = "TECCE",
    "gsynth" = "GSC",
    "gsynth with Control" = "GSC + $x_{it}$",
    "TWFE via OLS" = "TWFE",
    "TWFE with Control" = "TWFE + $x_{it}$"
  )
  res <- res |>
    mutate(
      estimator = table_order[
        match(estimator, names(table_order))
      ]
    ) |>
    arrange(
      match(estimator, table_order),
      N,
      T
    )

  res_wide <- res |>
    select(T, estimator, Bias = bias, RMSE = rmse) |>
    arrange(T) |>
    pivot_wider(
      id_cols = c("estimator"),
      names_from = "T",
      values_from = c(Bias, RMSE),
      names_glue = "{.value} - $T = {T}$",
      names_vary = "slowest"
    ) |>
    setNames(
      c("Estimator", rep(c("Bias", "RMSE"), 4))
    )

  tab <- res_wide |>
    tt(
      caption = tbl_caption
    ) |>
    format_tt(j = 2:ncol(res_wide), sprintf = "%0.3f") |>
    group_tt(
      j = list(
        "$T = 8$" = 2:3,
        "$T = 16$" = 4:5,
        "$T = 24$" = 6:7,
        "$T = 32$" = 8:9
      )
    )

  return(tab)
}

# %%
make_panels <- function(sim_results) {
  nested <- sim_results |>
    arrange(dgp_num) |>
    nest(
      .by = c(dgp_num, dgp_str, fact_label, eta_y, tau_x, avg_diff_in_loading)
    )

  tabs <- seq_len(nrow(nested)) |>
    map(function(i) {
      tbl_caption <- sprintf(
        "Factor specification: %s. DGP: %s.",
        nested$fact_label[i],
        nested$dgp_str[i]
      )

      nested$data[[i]] |>
        make_table(tbl_caption = tbl_caption)
    })

  fs::dir_create(here("out/tables/simulation"))
  seq_len(nrow(nested)) |>
    walk(function(i) {
      body_str <- tinytable:::build_tt(tabs[[i]], output = "latex")@body

      # 1, 2, and 3
      dgp_number <- 1 + ((nested$dgp_num[i] - 1) %% 3)
      fact_str <- janitor::make_clean_names(nested$fact_label[i])

      out_file <- here(
        "out/tables/simulation",
        sprintf("tab-long-fact_%s-dgp_%s.tex", fact_str, dgp_number)
      )
      cat(
        body_str,
        file = out_file
      )
    })
}


# %%
make_panels(sim_results)
