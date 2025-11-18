# %%
library(tidyverse)
library(tinytable)
library(here)
library(kfbmisc) # remotes::install_github("kylebutts/kfbmisc")

# fs::dir_create(here("out/tables/simulation/"))

# sim_results <- here("data/simulations/sim_results.csv") |>
#   read_csv(show_col_types = FALSE)

sim_results <- here("data/simulations-new/sim_results.csv") |>
  read_csv(show_col_types = FALSE)

sim_results <- sim_results |>
  mutate(
    fact_label = fct(
      fact_label,
      levels = c("Constant", "Linear Trend", "Autocorrelated")
    )
  ) |>
  mutate(
    abs_bias = abs(bias)
  )

# %%
plot_dgp_sim_results <- function(
  fact_label = c("Constant", "Linear Trend", "Autocorrelated"),
  avg_diff_in_loading = 0,
  tau_x = 0,
  eta_y = 0
) {
  subset <- sim_results |>
    filter(
      avg_diff_in_loading == .env$avg_diff_in_loading,
      fact_label == .env$fact_label,
      eta_y == .env$eta_y,
      tau_x == .env$tau_x
    ) |>
    arrange(T) |>
    mutate(
      position = case_when(
        estimator == "gsynth" ~ T - 0.2,
        estimator == "gsynth + unit FE" ~ T - 0.1,
        estimator == "TECCE" ~ T,
        estimator == "TWFE via OLS" ~ T + 0.1,
        estimator == "TWFE with Controls" ~ T + 0.2,
        TRUE ~ T
      )
    )

  ggplot(data = subset) +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = kfbmisc::tailwind_color("zinc-500")
    ) +
    geom_point(
      aes(x = position, y = bias, color = estimator),
    ) +
    geom_line(
      aes(x = position, y = bias, color = estimator, group = estimator),
    ) +
    scale_x_continuous(
      breaks = unique(sim_results$T)
    ) +
    scale_color_manual(
      values = c(
        kyle_color("navy"),
        kyle_color("purple"),
        kyle_color("green"),
        kyle_color("rose"),
        kyle_color("yellow"),
        kyle_color("magenta")
      ),
      guide = guide_legend(nrow = 1)
    ) +
    labs(
      x = "T",
      y = "Bias of estimator",
      color = NULL
    ) +
    kfbmisc::theme_kyle(base_size = 14, grid_minor = "h") +
    theme(
      legend.position = "bottom",
      legend.title.position = "top",
      legend.key.spacing.x = unit(12, "pt")
    )
}


# %%
plot_dgp_sim_results(
  fact_label = "Constant",
  avg_diff_in_loading = 0,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Constant",
  avg_diff_in_loading = -0.5,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Constant",
  avg_diff_in_loading = -0.5,
  tau_x = 1,
  eta_y = 0
)

# %%
plot_dgp_sim_results(
  fact_label = "Linear Trend",
  avg_diff_in_loading = 0,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Linear Trend",
  avg_diff_in_loading = -0.5,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Linear Trend",
  avg_diff_in_loading = -0.5,
  tau_x = 1,
  eta_y = 0
)

# %%
plot_dgp_sim_results(
  fact_label = "Autocorrelated",
  avg_diff_in_loading = 0,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Autocorrelated",
  avg_diff_in_loading = -0.5,
  tau_x = 0,
  eta_y = 0
)
plot_dgp_sim_results(
  fact_label = "Autocorrelated",
  avg_diff_in_loading = -0.5,
  tau_x = 1,
  eta_y = 0
)
