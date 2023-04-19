library(data.table)
library(here)
library(did2s)
library(ggplot2)
library(tikzDevice)

source(here("figures/convert_tex_to_pdf.R"))

# Load and clean data 

df <- haven::read_dta(here("code/Trade-Liberalization-and-Markup-Dispersion/data/AEJ_ind_DID_3-digit.dta"))
setDT(df)

# construct key variables
df[year == 2001, a01 := avgtariff_ind3 / 100]
df[, tariff01 := mean(a01[year == 2001]), by = "sic3"]
df[, post02 := (year > 2001)]
df[, t01post02 := tariff01 * post02]

median_tariff01 <- df[year == 2001, median(tariff01, na.rm = TRUE)]
df[,
  hightariff := (tariff01[year == 2001] >= median_tariff01),
  by = "sic3"
]
df[, hightariff_post02 := hightariff & post02]
df[, rel_year := ifelse(hightariff, year - 2002, -Inf)]

df[, ln_gini := log(gini)]
df[, ln_theil := log(theil)]
df[, ln_cv := log(cv)]
df[, ln_mld := log(mld)]
df[, ln_rmd := log(rmd)]

df[, lnn := log(n)]
df[, lnasset := log(assets)]
df[, lnexports := log(exports)]
df[, lnfdi := log(foreign)]
df[, lnit := log(input_tariff)]
df[, tariff := avgtariff_ind3 / 100]

df <- df[!is.na(hightariff), ]

# N = 164
# df$sic3 |> unique() |> length()
# T = 9
# df$year |> unique() |> length()


# Figure 4 ---------------------------------------------------------------------

# feols(
#   theil ~ 0 + i(year, ref = 1998) | sic3,
#   df,
#   split = ~hightariff
# ) |>
#   iplot()


# TWFE -------------------------------------------------------------------------

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

# iplot(
#   list(est_ln_theil_twfe, est_ln_theil_twfe_cov)
# )

# did2s ------------------------------------------------------------------------

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

# coefplot(list(est_ln_theil_did2s, est_ln_theil_did2s_cov))



# CCE DID ----------------------------------------------------------------------

setorder(df, sic3, year)
setDT(df)

df$y <- df$ln_theil
df$x1 <- df$eg_3dt_city
df$x2 <- df$ln_theil_tfp
# df$x3 <- 1
# str(df)

xvars <- c("x1", "x2")

df <- df[!(is.na(y) | is.na(x1) | is.na(x2)), ]
df[, n := .N, by = "sic3"]
df <- df[n == 8, ]
T0 <- 2001
df[, g := ifelse(hightariff == TRUE, 2002, Inf)]

# This is a lazy way of not doing the within-transformation without needing to modify all the code
df$ytilde <- df$y
df$x1tilde <- df$x1
df$x2tilde <- df$x2
# df$x3tilde <- df$x3
# df$x4tilde <- df$x4

# Export data
# fwrite(
#   df,
#   file = here("code/Trade-Liberalization-and-Markup-Dispersion/data/sample.csv")
# )

## CCE pooled estimate of \beta hat --------------------------------------------

Fhat <- df |>
  DT(
    g == Inf,
    lapply(.SD, \(x) mean(x)),
    by = year,
    .SDcols = c("ytilde", paste0(xvars, "tilde"))
  ) |>
  DT(,
    as.matrix(.SD),
    .SDcols = c("ytilde", paste0(xvars, "tilde"))
  )

# Add sic3-specific linear-time trends
# Fhat = cbind(Fhat, 1:nrow(Fhat) * 1e-12)


N_pre_T0 <- (T0 - min(df$year) + 1)
Fpre <- Fhat[1:N_pre_T0, ]
# M_Fpre = diag(N_pre_T0) - Fpre %*% solve(t(Fpre) %*% Fpre) %*% t(Fpre)
M_Fpre <- diag(N_pre_T0) - Fpre %*% solve(t(Fpre) %*% Fpre) %*% t(Fpre)

B <- matrix(0, nrow = length(xvars), ncol = length(xvars))
A <- matrix(0, nrow = length(xvars), ncol = 1)

for (id in unique(df[, sic3])) {
  Xi_pre <- df |>
    DT(
      sic3 == id & year <= T0,
      as.matrix(.SD),
      .SDcols = paste0(xvars, "tilde")
    )

  yi_pre <- df[
    sic3 == id & year <= T0,
    ytilde
  ]

  B <- B + t(Xi_pre) %*% M_Fpre %*% Xi_pre
  A <- A + t(Xi_pre) %*% M_Fpre %*% yi_pre
}

# CCEP estimator using pre-treatment X's for all groups
bhat <- solve(B) %*% A


## Imputation of covariates and outcome ----------------------------------------

for (id in unique(df[g < Inf, sic3])) {
  # impute factor loadings
  Xi_pre <- df |>
    DT(
      sic3 == id & year <= T0,
      as.matrix(.SD),
      .SDcols = paste0(xvars, "tilde")
    )

  yi_pre <- df[
    sic3 == id & year <= T0,
    ytilde
  ]

  ghats <- solve(t(Fpre) %*% Fpre) %*% t(Fpre) %*% (yi_pre - Xi_pre %*% bhat)

  # impute X(0)
  X_0hat <- Fhat %*% solve(t(Fpre) %*% Fpre) %*% t(Fpre) %*% Xi_pre

  for (i in 1:length(xvars)) {
    var_name <- paste0(xvars[i], "tilde_0hat")
    df |>
      DT(
        sic3 == id,
        var_name := X_0hat[, i],
        env = list(var_name = var_name)
      )
  }

  # impute y(0) using X(0)
  y_0hat <- (X_0hat %*% bhat) + (Fhat %*% ghats)

  df[sic3 == id, ytilde_0hat := y_0hat]

  # impute y(0) using observed X
  Xi <- df |>
    DT(
      sic3 == id,
      as.matrix(.SD),
      .SDcols = paste0(xvars, "tilde")
    )

  y_0hat_obs_x <- (Xi %*% bhat) + (Fhat %*% ghats)

  df[sic3 == id, ytilde_0hat_obs_x := y_0hat_obs_x]
}

# Difference
df$y_diff <- df$ytilde - df$ytilde_0hat
df$y_diff_obs_x <- df$ytilde - df$ytilde_0hat_obs_x
df$x1_diff <- (df$x1tilde - df$x1tilde_0hat)
df$mediated1_diff <- df$x1_diff * bhat[1]
df$x2_diff <- (df$x2tilde - df$x2tilde_0hat)
df$mediated2_diff <- df$x2_diff * bhat[2]


## Estimate Effects ------------------------------------------------------------

# Log(theil)
est_ln_theil_cce <- feols(
  y_diff ~ 0 + i(rel_year),
  df[g < Inf, ]
)

est_ln_theil_cce_obs_x <- feols(
  y_diff_obs_x ~ 0 + i(rel_year),
  df[g < Inf, ]
)

# [X_{1it} - X_{1it}(0)] * \beta_1
est_x1_cce <- feols(
  mediated1_diff ~ 0 + i(rel_year),
  df[g < Inf, ]
)

# [X_{2it} - X_{2it}(0)] * \beta_2
est_x2_cce <- feols(
  mediated2_diff ~ 0 + i(rel_year),
  df[g < Inf, ]
)


# Plot Results ----------------------------------------------------------------

coefplot(
  list(est_ln_theil_did2s, est_ln_theil_cce, est_ln_theil_cce_obs_x),
  keep = "rel_year"
)
# coefplot(
#   list(est_ln_theil_did2s, est_ln_theil_cce, est_x1_cce, est_x2_cce),
#   keep = "rel_year"
# )
# coefplot(est_x1_cce)
# coefplot(est_x2_cce)

## TWFE Imputation vs. CCEDID -------------------------------------------------

est1 <- broom::tidy(est_ln_theil_did2s)
est1$estimator <- "TWFE Imputation"
est1$estimator_num <- 1
est2 <- broom::tidy(est_ln_theil_cce)
est2$estimator <- "CCEDID"
est2$estimator_num <- 2
# est3 <- broom::tidy(est_x2_cce)
# est3$estimator <- "CCEDID Mechanism"
# est3$estimator_num <- 3

ests <- rbindlist(list(est1, est2))
setnames(ests, "std.error", "se")
ests$rel_year <- stringr::str_replace(ests$term, "rel_year::", "") |>
  as.numeric()

ests[, rel_year := data.table::fcase(
  estimator_num == 1, rel_year - 0.1,
  estimator_num == 2, rel_year + 0.1
)]

ests$estimator <- forcats::as_factor(ests$estimator)

# https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")

(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.18, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se,
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
    size = 2.5
  ) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  scale_color_manual(values = c(colors[7], "black")) +
  scale_shape_manual(values = c(15, 16, 18, 15)) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)', color = NULL, shape = NULL) +
  guides(color = guide_legend(
    # label.position = "bottom",
    override.aes = list(linetype = 0)
  )) +
  theme_bw(base_size = 9) +
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


# ggsave(here("figures/trade-est.pdf"), p, width = 10 * 6.9/10, height = 5 * 6.9/10, bg = "white")

tikzDevice::tikz(
  here("figures/trade-est.tex"), 
  width = 10 * 6.9 / 10, height = 5 * 6.9 / 10,
  standAlone = FALSE
)
plot(p)
dev.off()
compile_tikz(here("figures/trade-est.tex"))

## Decomposition of Effects ---------------------------------------------------

est1 <- broom::tidy(est_ln_theil_cce)
est1$estimator <- r'(Total effect: $\Delta_t + \gamma_{1t} + \gamma_{2t}$)'
est1$estimator_num <- 1

est2 <- broom::tidy(est_ln_theil_cce_obs_x)
est2$estimator <- r'(Direct effect: $\Delta_t$)'
est2$estimator_num <- 2

est3 <- broom::tidy(est_x1_cce)
est3$estimator <- r'(Indirect effect: $\gamma_{1t}$)'
est3$estimator_num <- 3

est4 <- broom::tidy(est_x2_cce)
est4$estimator <- r'(Indirect effect: $\gamma_{2t}$)'
est4$estimator_num <- 4

ests <- rbindlist(list(est1, est2, est3, est4))
setnames(ests, "std.error", "se")
ests$rel_year <- stringr::str_replace(ests$term, "rel_year::", "") |>
  as.numeric()

ests[, rel_year := data.table::fcase(
  estimator_num == 1, rel_year - 0.15,
  estimator_num == 2, rel_year + 0.1,
  estimator_num == 3, rel_year + 0.2,
  estimator_num == 4, rel_year + 0.3
)]

ests$estimator <- forcats::as_factor(ests$estimator)

# https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")



(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.2, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se,
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
    size = 2
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
  theme_bw(base_size = 9) +
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

# ggsave(here("figures/trade-est.pdf"), p, width = 10 * 6.9/10, height = 5 * 6.9/10, bg = "white")

tikzDevice::tikz(here("figures/trade-decomposition.tex"), width = 10 * 6.9 / 10, height = 5 * 6.9 / 10, bg = "white", standAlone = FALSE)
plot(p)
dev.off()
compile_tikz(here("figures/trade-decomposition.tex"))



## Using observed X vs. X(0) ---------------------------------------------------

est1 <- broom::tidy(est_ln_theil_cce)
est1$estimator <- r'(CCEDID using $\hat{X}_{it}(0)$)'
est1$estimator_num <- 1
est2 <- broom::tidy(est_ln_theil_cce_obs_x)
est2$estimator <- r'(CCEDID using $X_{it}$)'
est2$estimator_num <- 2
est3 <- broom::tidy(est_x2_cce)
est3$estimator <- "CCEDID Mechanism"
est3$estimator_num <- 3

ests <- rbindlist(list(est1, est2, est3))
setnames(ests, "std.error", "se")
ests$rel_year <- stringr::str_replace(ests$term, "rel_year::", "") |>
  as.numeric()

ests[, rel_year := data.table::fcase(
  estimator_num == 1, rel_year - 0.15,
  estimator_num == 2, rel_year,
  estimator_num == 3, rel_year + 0.15
)]

ests$estimator <- forcats::as_factor(ests$estimator)

# https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")

(p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2002, linetype = "dashed") +
  annotate(
    "text",
    x = 2002.1, y = 0.2, hjust = 0,
    label = r'($\leftarrow$ China joins WTO)', size = 3
  ) +
  geom_errorbar(
    data = ests,
    aes(
      x = rel_year + 2002,
      ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se,
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
    size = 2.5
  ) +
  scale_x_continuous(breaks = 2002 + -4:3) +
  scale_color_manual(values = c(colors[4], colors[8], colors[1], colors[5])) +
  scale_shape_manual(values = c(15, 16, 18, 15)) +
  labs(x = NULL, y = r'(Estimated Effect on $\log($Theil$)$)', color = NULL, shape = NULL) +
  guides(color = guide_legend(
    # label.position = "bottom",
    override.aes = list(linetype = 0)
  )) +
  theme_bw(base_size = 9) +
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

tikzDevice::tikz(here("figures/trade-x0-vs-obs-x.tex"), width = 10 * 6.9 / 10, height = 5 * 6.9 / 10, bg = "white", standAlone = FALSE)
plot(p)
dev.off()
compile_tikz(here("figures/trade-x0-vs-obs-x.tex"))


