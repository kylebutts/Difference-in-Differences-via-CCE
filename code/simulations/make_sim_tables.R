# %%
library(tidyverse)
library(tinytable)
library(here)
library(kfbmisc) # remotes::install_github("kylebutts/kfbmisc")

# sim_results <- here("data/simulations/sim_results.csv") |>
#   read_csv(show_col_types = FALSE)

make_table <- function(res) {
  estimator_order <- c(
    "TECCE",
    "gsynth",
    "gsynth with Control",
    "gsynth with TWFE",
    "TWFE via OLS",
    "TWFE with Controls"
  )
  res <- res |>
    arrange(
      match(estimator, estimator_order),
      N,
      T
    )

  res_wide <- res |>
    select(N, T, estimator, Bias = bias, RMSE = rmse) |>
    pivot_wider(
      id_cols = c("N", "T"),
      names_from = estimator,
      values_from = c(Bias, RMSE),
      names_sep = " - "
    )

  tab <- res_wide |>
    tt() |>
    format_tt(j = 3:ncol(res_wide), sprintf = "%0.3f")

  tab |>
    tinytable:::build_tt("latex") |>
    (\(x) x@body)()
}

# %%
make_panels <- function(sim_results) {
  results_wider <- sim_results |>
    select(
      N,
      T,
      dgp_number,
      dgp_label,
      fact_label,
      err_label,
      eta_y,
      tau_x,
      avg_diff_in_loading,
      estimator,
      Bias = bias,
      RMSE = rmse
    ) |>
    pivot_longer(
      cols = c(Bias, RMSE),
      names_to = "statistic",
      values_to = "value"
    ) |>
    pivot_wider(
      id_cols = c(
        N,
        T,
        dgp_number,
        dgp_label,
        fact_label,
        err_label,
        eta_y,
        tau_x,
        avg_diff_in_loading,
        statistic
      ),
      names_from = estimator,
      values_from = c(value)
    )

  nested <- results_wider |>
    arrange(
      match(fact_label, c("Constant", "Linear Trend", "Autocorrelated"))
    ) |>
    nest(.by = c(dgp_number, dgp_label, fact_label, err_label))

  tabs <- seq_len(nrow(nested)) |>
    map(function(i) {
      tbl_caption <- sprintf(
        "Factor specification: %s. DGP: %s. Error term: %s.",
        nested$fact_label[i],
        nested$dgp_label[i],
        nested$err_label[i]
      )

      tab_data <- nested$data[[i]]
      tab_data <- tab_data |>
        filter(N != 100)
      rmse_rows <- tab_data$statistic == "RMSE"
      tab_data <- tab_data |>
        arrange(N, T) |>
        select(-eta_y, -tau_x, -avg_diff_in_loading)

      tab_data <- tab_data |> select(-statistic)
      tab_data |>
        tt(
          caption = tbl_caption
        ) |>
        format_tt(
          i = which(!rmse_rows),
          j = 3:ncol(tab_data),
          sprintf = "%0.2f"
        ) |>
        format_tt(
          i = which(rmse_rows),
          j = 3:ncol(tab_data),
          sprintf = "[%0.3f]"
        ) |>
        format_tt(
          i = which(rmse_rows),
          j = c("N", "T"),
          fn = function(x) {
            ""
          }
        )
    })

  ## Export
  fs::dir_create(here("out/tables/simulation"))
  seq_len(nrow(nested)) |>
    walk(function(i) {
      body_str <- tinytable:::build_tt(tabs[[i]], output = "latex")@body

      # 1, 2, and 3
      dgp_number <- nested$dgp_number[i]
      fact_str <- janitor::make_clean_names(nested$fact_label[i])
      is_autocorrelated_error <- +(nested$err_label[i] == "AR(1)")

      out_file <- here(
        "out/tables/simulation",
        sprintf(
          "tab-fact_%s-dgp_%s-eps_auto_%s.tex",
          fact_str,
          dgp_number,
          is_autocorrelated_error
        )
      )
      cat(
        body_str,
        file = out_file
      )
    })
}


# %%
make_panels(sim_results)
