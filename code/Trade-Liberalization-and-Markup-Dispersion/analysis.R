#' ---
#' title: "Main analysis"
#' ---

# %%
#| warning: false
library(data.table)
library(tidyverse)
library(here)
library(did2s)
library(ggplot2)
library(kfbmisc)
library(fixest)
library(stringmagic)
setFixest_etable(markdown = TRUE)

# https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")

# N bootstraps
R <- 499

# Load and clean data ----------------------------------------------------------
# %%
df <- haven::read_dta(here("code/Trade-Liberalization-and-Markup-Dispersion/data/AEJ_ind_DID_3-digit.dta"))

# construct key variables
df <- df |>
  mutate(
    tariff01 = ifelse(
      any(year == 2001),
      avgtariff_ind3[year == 2001] / 100, NA_real_
    ),
    .by = "sic3"
  ) |>
  mutate(
    hightariff = tariff01 > median(tariff01, na.rm = TRUE),
    g = ifelse(hightariff == TRUE, 2002, Inf),
    hightariff_post02 = hightariff & year > 2001,
    rel_year = if_else(hightariff, year - 2002, -Inf),
    pre_post = case_when(
      rel_year > -Inf & rel_year < 0 ~ "pre",
      rel_year >= 0 ~ "post",
      .default = "control"
    )
  ) |>
  filter(!is.na(tariff01))

df <- df |>
  mutate(
    ln_gini = log(gini),
    ln_theil = log(theil),
    ln_cv = log(cv),
    ln_mld = log(mld),
    ln_rmd = log(rmd),
    lnn = log(n),
    lnasset = log(assets),
    lnexports = log(exports),
    lnfdi = log(foreign),
    lnit = log(input_tariff),
    tariff = avgtariff_ind3 / 100
  )

# Figure 4 ---------------------------------------------------------------------
# %%
feols(
  theil ~ 0 + i(year, ref = 1998) | sic3,
  df,
  split = ~hightariff
) |>
  iplot()

# TWFE -------------------------------------------------------------------------
# %%
est_ln_theil_twfe <- feols(
  ln_theil ~ i(rel_year, ref = c(-1, -Inf)) | sic3 + year,
  df,
  cluster = ~sic3
)

# Add time-varying covariates
est_ln_theil_twfe_cov <- feols(
  ln_theil ~ ln_theil_tfp + eg_3dt_city + i(rel_year, ref = c(-1, -Inf)) | sic3 + year,
  df,
  cluster = ~sic3
)

# did2s ------------------------------------------------------------------------
# %%
est_ln_theil_did2s <- did2s:::did2s(
  df,
  yname = "ln_theil",
  first_stage = ~ 0 | sic3 + year,
  second_stage = ~ i(rel_year, ref = -Inf),
  treatment = "hightariff_post02",
  cluster_var = "sic3"
)

# Add time-varying covariates
est_ln_theil_did2s_cov <- did2s:::did2s(
  df,
  yname = "ln_theil",
  first_stage = ~ ln_theil_tfp + eg_3dt_city | sic3 + year,
  second_stage = ~ i(rel_year, ref = -Inf),
  treatment = "hightariff_post02",
  cluster_var = "sic3"
)

# %%
est_ln_theil_did2s_pre_post <- did2s:::did2s(
  df,
  yname = "ln_theil",
  first_stage = ~ ln_theil_tfp + eg_3dt_city | sic3 + year,
  second_stage = ~ i(pre_post, ref = "control"),
  treatment = "hightariff_post02",
  cluster_var = "sic3"
)


# CCE DID ----------------------------------------------------------------------
# %%
df <- df |>
  arrange(sic3, year) |>
  mutate(y = ln_theil, x1 = ln_theil_tfp) |>
  filter(!is.na(y), !is.na(x1)) |>
  # Balanced panel
  filter(n() == 8, .by = "sic3")

T0 <- 2001
xvars <- c("x1")

# %%
collapse <- df |>
  summarize(
    mean_ln_theil_tfp = mean(ln_theil_tfp, na.rm = TRUE),
    sd_ln_theil_tfp = sd(ln_theil_tfp, na.rm = TRUE),
    .by = year
  )

plot_ln_theil_tfp <- ggplot() +
  geom_point(
    aes(x = year, y = ln_theil_tfp),
    data = df,
    size = 1.1
  ) +
  geom_point(
    aes(x = year, y = mean_ln_theil_tfp),
    data = collapse,
    color = "blue", size = 3, shape = 16
  ) +
  labs(x = NULL, y = r'($\ln($Theil TFP$)$)') +
  kfbmisc::theme_kyle(base_size = 10)


## CCE pooled estimate of \beta hat --------------------------------------------
# %%
# CSAs of y and Xs
Fhat <- df |>
  filter(g == Inf) |>
  summarize(
    across(c(y, all_of(xvars)), ~ mean(.x, na.rm = TRUE)),
    .by = "year"
  ) |>
  select(-year) |>
  as.matrix() |>
  # Add constant to Fhat (unit fixed-effects)
  cbind(1)

stringmagic::cat_magic("The estimated column rank is {qr(Fhat)$rank} out of {ncol(Fhat)} columns")

# CCEP estimator for \hat{\beta} using pre-treatment X's for all groups
N_pre_T0 <- (T0 - min(df$year) + 1)
Fpre <- Fhat[1:N_pre_T0, ]
M_Fpre <- diag(N_pre_T0) - Fpre %*% solve(crossprod(Fpre, Fpre), t(Fpre))

B <- matrix(0, nrow = length(xvars), ncol = length(xvars))
A <- matrix(0, nrow = length(xvars), ncol = 1)
for (id in unique(df$sic3)) {
  idx <- (df$sic3 == id) & (df$year <= T0)
  Xi_pre <- as.matrix(df[idx, xvars])
  yi_pre <- df[idx, ]$y

  B <- B + t(Xi_pre) %*% M_Fpre %*% Xi_pre
  A <- A + t(Xi_pre) %*% M_Fpre %*% yi_pre
}
bhat <- solve(B, A)

# %%
impute_id <- function(sub) {
  # impute factor loadings
  Xi_pre <- sub |>
    filter(year <= T0) |>
    select(all_of(xvars)) |>
    as.matrix()

  Xi <- sub |>
    select(all_of(xvars)) |>
    as.matrix()

  yi_pre <- sub |>
    filter(year <= T0) |>
    select(y) |>
    as.matrix()

  # For total effect, do not need to impute X
  ghats_ignore_X <- solve(crossprod(Fpre), crossprod(Fpre, yi_pre))
  y_0hat_ignore_X <- (Fhat %*% ghats_ignore_X)

  # Estimate factor-loadings
  ghats <- solve(crossprod(Fpre), crossprod(Fpre, yi_pre - Xi_pre %*% bhat))

  # impute X(0)
  X_0hat <- Fhat %*% solve(crossprod(Fpre), crossprod(Fpre, Xi_pre))
  colnames(X_0hat) <- paste0(xvars, "_0hat")

  # impute y(0) using X(0)
  y_0hat <- (X_0hat %*% bhat) + (Fhat %*% ghats)

  # impute y(0) using observed X
  y_0hat_obs_x <- (Xi %*% bhat) + (Fhat %*% ghats)

  imputed <- tibble(
    y_0hat_ignore_X = as.numeric(y_0hat_ignore_X),
    x1_0hat = as.numeric(X_0hat),
    y_0hat = as.numeric(y_0hat),
    y_0hat_obs_x = as.numeric(y_0hat_obs_x)
  )
  return(imputed)
}

df <- df |>
  mutate(
    impute_id(pick(everything())),
    .by = "sic3"
  )

df <- df |>
  mutate(
    y_diff = y - y_0hat,
    y_diff_ignore_X = y - y_0hat_ignore_X,
    y_diff_obs_x = y - y_0hat_obs_x,
    x1_diff = (x1 - x1_0hat),
    mediated1_diff = x1_diff * bhat[1]
  )

## Estimate Effects ------------------------------------------------------------

# %%
# Showing equivalence
etable(feols(
  c(y_diff, y_diff_ignore_X) ~ 0 + i(rel_year),
  df |> filter(g != Inf),
  vcov = "hc1"
))

# %%
# Log(theil)
est_ln_theil_cce <- feols(
  y_diff ~ 0 + i(rel_year),
  df |> filter(g != Inf),
  vcov = "hc1"
)
res_bootstrap <- boot::boot(
  df |> filter(g != Inf),
  function(data, indices) {
    coef(feols(
      y_diff ~ 0 + i(rel_year),
      data[indices, ]
    ))
  },
  R = R,
  strata = df |> filter(g != Inf) |> pull(sic3)
)
est_ln_theil_cce <- summary(est_ln_theil_cce, cov(res_bootstrap$t))

# %%
# [X_{1it} - X_{1it}(0)] * \beta_1
est_x1_cce <- feols(
  mediated1_diff ~ 0 + i(rel_year),
  df |> filter(g != Inf),
  vcov = "hc1"
)
res_bootstrap <- boot::boot(
  df |> filter(g != Inf),
  function(data, indices) {
    coef(feols(
      mediated1_diff ~ 0 + i(rel_year),
      data[indices, ]
    ))
  },
  R = R,
  strata = df |> filter(g != Inf) |> pull(sic3)
)
est_x1_cce <- summary(est_x1_cce, cov(res_bootstrap$t))

# %%
est_ln_theil_cce_obs_x <- feols(
  y_diff_obs_x ~ 0 + i(rel_year),
  df |> filter(g != Inf),
  vcov = "hc1"
)
res_bootstrap <- boot::boot(
  df |> filter(g != Inf),
  function(data, indices) {
    coef(feols(
      y_diff_obs_x ~ 0 + i(rel_year),
      data[indices, ]
    ))
  },
  R = R,
  strata = df |> filter(g != Inf) |> pull(sic3)
)
est_ln_theil_cce_obs_x <- summary(est_ln_theil_cce_obs_x, cov(res_bootstrap$t))

# %%
est_ln_theil_pre_post <- feols(
  y_diff ~ 0 + i(pre_post),
  df |> filter(g != Inf),
  vcov = "hc1"
)
est_x1_pre_post <- feols(
  mediated1_diff ~ 0 + i(pre_post),
  df |> filter(g != Inf),
  vcov = "hc1"
)

# Plot Results ----------------------------------------------------------------
## Graph: CCEDID ---------------------------------------------------------------
# %%
ests <- broom::tidy(est_ln_theil_cce) |>
  mutate(
    rel_year = as.numeric(str_replace(term, "rel_year::", ""))
  )

ests_pre_post <- broom::tidy(est_ln_theil_pre_post) |>
  mutate(pre_post = str_replace(term, "pre_post::", ""))

ests_pre_post <- bind_rows(
  ests_pre_post |> filter(pre_post == "pre"),
  ests_pre_post |> filter(pre_post == "pre"),
  ests_pre_post |> filter(pre_post == "post"),
  ests_pre_post |> filter(pre_post == "post")
)
ests_pre_post$x <- c(-4, -1, 0, 3)

# %%
(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.15, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error
    ),
    color = colors[1],
    # Change to 1 for pdf version
    linewidth = 1,
    width = 0.08
  ) +
  geom_point(
    data = ests,
    aes(
      x = rel_year + 2002, y = estimate
    ),
    size = 1,
    color = colors[1], shape = 15
  ) +
  geom_ribbon(
    data = ests_pre_post,
    aes(
      x = x + 2002,
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error,
      group = pre_post
    ),
    fill = colors[1],
    alpha = 0.1,
    color = colors[1],
    linetype = "dotted",
  ) +
  geom_line(
    data = ests_pre_post,
    aes(
      x = x + 2002,
      y = estimate,
      group = pre_post
    ),
    color = colors[1],
    linewidth = 1,
    linetype = "dotted"
  ) +
  scale_y_continuous(limits = c(-0.5, 0.15)) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)') +
  kfbmisc::theme_kyle(base_size = 10) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = -10)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/trade-cce_est.pdf"),
  p,
  width = 7, height = 3.5
)

## Graph: X1 Mediated ----------------------------------------------------------
# %%
ests <- broom::tidy(est_x1_cce) |>
  mutate(
    rel_year = as.numeric(str_replace(term, "rel_year::", ""))
  )

ests_pre_post <- broom::tidy(est_x1_pre_post) |>
  mutate(pre_post = str_replace(term, "pre_post::", ""))

ests_pre_post <- bind_rows(
  ests_pre_post |> filter(pre_post == "pre"),
  ests_pre_post |> filter(pre_post == "pre"),
  ests_pre_post |> filter(pre_post == "post"),
  ests_pre_post |> filter(pre_post == "post")
)
ests_pre_post$x <- c(-4, -1, 0, 3)

# %%
(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.15, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  annotate(
    geom = "label", x = 1998, y = -0.07,
    label = r'($\hat{\beta} = 0.238$)',
    bg = "white", hjust = 0,
    label.padding = unit(0.5, "lines"),
    label.r = unit(0, "lines"),
    label.size = 0
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error
    ),
    color = colors[7],
    # Change to 1 for pdf version
    linewidth = 1,
    width = 0.08
  ) +
  geom_point(
    data = ests,
    aes(
      x = rel_year + 2002, y = estimate
    ),
    size = 1,
    color = colors[7], shape = 15
  ) +
  geom_ribbon(
    data = ests_pre_post,
    aes(
      x = x + 2002,
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error,
      group = pre_post
    ),
    fill = colors[7],
    alpha = 0.1,
    color = colors[7],
    linetype = "dotted",
  ) +
  geom_line(
    data = ests_pre_post,
    aes(
      x = x + 2002,
      y = estimate,
      group = pre_post
    ),
    color = colors[7],
    linewidth = 1,
    linetype = "dotted"
  ) +
  scale_y_continuous(limits = c(-0.1, 0.15)) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)') +
  kfbmisc::theme_kyle(base_size = 10) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = -10)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/trade-cce_mediated_est.pdf"),
  p,
  width = 7, height = 3.5
)

## Graph: Using observed X vs. X(0) --------------------------------------------
# %%
est1 <- broom::tidy(est_ln_theil_cce)
est1$estimator <- r'(CCEDID using $\hat{X}_{it}(0)$)'
est1$estimator_num <- 1
est2 <- broom::tidy(est_ln_theil_cce_obs_x)
est2$estimator <- r'(CCEDID using $X_{it}$)'
est2$estimator_num <- 2

ests <- bind_rows(est1, est2) |>
  mutate(
    rel_year = str_replace(term, "rel_year::", ""),
    rel_year = as.numeric(rel_year),
    rel_year = case_when(
      estimator_num == 1 ~ rel_year - 0.075,
      estimator_num == 2 ~ rel_year + 0.075
    ),
    estimator = as_factor(estimator)
  )

# %%
(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.15, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error,
      color = estimator
    ),
    # Change to 1 for pdf version
    linewidth = 1.5,
    width = 0.08
  ) +
  geom_point(
    data = ests,
    aes(
      x = rel_year + 2002, y = estimate,
      color = estimator, shape = estimator
    ),
    size = 1.5
  ) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  scale_color_manual(values = c(colors[1], colors[2])) +
  scale_shape_manual(values = c(15, 16)) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)', color = NULL, shape = NULL) +
  guides(color = guide_legend(
    # label.position = "bottom",
    override.aes = list(linetype = 0)
  )) +
  theme_bw(base_size = 10) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = -10)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/trade-x0-vs-obs-x.pdf"),
  p,
  width = 7, height = 3.5
)

## Graph: TWFE Imputation vs. CCEDID -------------------------------------------
# %%
est1 <- broom::tidy(est_ln_theil_did2s)
est1$estimator <- "TWFE Imputation"
est1$estimator_num <- 1
est2 <- broom::tidy(est_ln_theil_cce)
est2$estimator <- "CCEDID"
est2$estimator_num <- 2

ests <- bind_rows(est1, est2) |>
  mutate(
    rel_year = str_replace(term, "rel_year::", ""),
    rel_year = as.numeric(rel_year),
    rel_year = case_when(
      estimator_num == 1 ~ rel_year - 0.075,
      estimator_num == 2 ~ rel_year + 0.075
    ),
    estimator = as_factor(estimator)
  )

# %%
(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.15, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error,
      color = estimator
    ),
    # Change to 1 for pdf version
    linewidth = 1.5,
    width = 0.08
  ) +
  geom_point(
    data = ests,
    aes(
      x = rel_year + 2002, y = estimate,
      color = estimator, shape = estimator
    ),
    size = 1.5
  ) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  scale_color_manual(values = c(colors[7], "black")) +
  scale_shape_manual(values = c(15, 16, 18, 15)) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)', color = NULL, shape = NULL) +
  guides(color = guide_legend(
    # label.position = "bottom",
    override.aes = list(linetype = 0)
  )) +
  theme_bw(base_size = 10) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = -10)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/trade-est.pdf"),
  p,
  width = 7, height = 3.5
)

## Graph: Decomposition of Effects ---------------------------------------------
# %%
est1 <- broom::tidy(est_ln_theil_cce)
est1$estimator <- r'(Total effect: $\Delta_t + \gamma_{1t}$)'
est1$estimator_num <- 1

est2 <- broom::tidy(est_ln_theil_cce_obs_x)
est2$estimator <- r'(Direct effect: $\Delta_t$)'
est2$estimator_num <- 2

est3 <- broom::tidy(est_x1_cce)
est3$estimator <- r'(Indirect effect: $\gamma_{1t}$)'
est3$estimator_num <- 3

ests <- bind_rows(est1, est2, est3) |>
  mutate(
    rel_year = str_replace(term, "rel_year::", ""),
    rel_year = as.numeric(rel_year),
    rel_year = case_when(
      estimator_num == 1 ~ rel_year - 0.15,
      estimator_num == 2 ~ rel_year + 0.1,
      estimator_num == 3 ~ rel_year + 0.2,
      estimator_num == 4 ~ rel_year + 0.3
    ),
    estimator = as_factor(estimator)
  )

# %%
(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.15, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error,
      color = estimator,
      alpha = estimator
    ),
    # Change to 1 for pdf version
    linewidth = 1.5,
    width = 0.08
  ) +
  geom_point(
    data = ests,
    aes(
      x = rel_year + 2002, y = estimate,
      color = estimator, shape = estimator
    ),
    size = 1.5
  ) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  scale_color_manual(values = c("black", "gray35", "gray35", "gray35")) +
  scale_shape_manual(values = c(15, 16, 18, 18)) +
  scale_alpha_manual(values = c(1, 1, 1, 1)) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)', color = NULL, shape = NULL, alpha = NULL) +
  guides(color = guide_legend(
    # label.position = "bottom",
    override.aes = list(linetype = 0, alpha = 1)
  )) +
  theme_bw(base_size = 10) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = -10)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/trade-decomposition.pdf"),
  p,
  width = 7, height = 3.5
)


# Cross-correlation tests ------------------------------------------------------
# Pesaran - 2004 - General Diagnostic Tests for Cross Section Dependence in Panels
source(here("code/Trade-Liberalization-and-Markup-Dispersion/utils.R"))

# "Defactorized"
eps <- df |>
  # filter(g != Inf) |>
  # identical with `y_diff`
  with(make_mat(year, sic3, y_diff_ignore_X, names = TRUE))

pesaran_test(eps)

# Manual: note this is not working well b/c eps is not mean zero (by column)
# eq (9)
# xi_it = e_{it} / || e_i ||_2
xi <- apply(eps, 2, \(x) x / norm(x, "2"))

# eq (8)
# rho_ij = \sum_t xi_it xi_jt
rho <- crossprod(xi)

# eq (7)
N <- ncol(xi)
T <- nrow(xi)
CD <- sqrt((2 * T) / (N * (N - 1))) * sum(rho[lower.tri(rho)])
p.value <- 2 * (1 - pnorm(abs(CD)))


# Testing rank condition (Da Vos et. al., 2014) --------------------------------
# %%
# Create Z
source(here("code/Trade-Liberalization-and-Markup-Dispersion/utils.R"))
temp <- df |> 
  arrange(sic3, year) |> 
  filter(hightariff == 0, year <= 2001)
vars <- c("y", "x1")
Z <- lapply(unique(temp$sic3), function(s) {
  as.matrix(temp[temp$sic3 == s, c("y", "x1")])
})

# Number of factors in CSAs
m_hat <- DVES_m(Z)
m_hat_y <- DVES_m(lapply(Z, function(Zi) Zi[, 1, drop = FALSE]))

# Rank of cross-sectional averages
rho_hat <- DVES_rho(Z)
cat_magic(
  "Estimated number of factors (m): {m_hat}\n", 
  "Estimated rank of CCE matrix (rho): {rho_hat}\n"
)

