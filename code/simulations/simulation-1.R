library(data.table)
library(ggplot2)
library(fixest)
library(did2s)

dgp <- 3
parallel_trends <- TRUE
# genate AR(1) error with rho = 0.75
set.seed(1)
T <- 8
# f2 <- as.numeric(arima.sim(model = list(ar = 0.8), n = T, n.start = 100))
# f2 <- 1:T
f2 <- c(1, 1, 1, 1, 1, 1, 1, 1)
params <- list(
  N = 200, T = T, T0 = 5,
  f2 = f2,
  b = c(1, 1), k = 2, 
  tau = 1, Delta = 1
)
S <- 2

simulate_data <- function(params, dgp, parallel_trends) {
  N <- params$N
  T <- params$T
  T0 <- params$T0
  b <- params$b
  K <- params$K

  f1 <- rep(1, T)
  f2 <- params$f2
  f <- cbind(f1, f2)
  p <- 2

  # Effect on y and X
  Delta <- params$Delta * (1:T > T0) * (dgp >= 2)
  tau <- params$tau * (1:T > T0) * (dgp == 3)

  # Correlation matrix for errors
  C <- matrix(0, nrow = T, ncol = T)
  for (t in 1:T) {
    for (s in 1:T) {
      C[t, s] <- 0.75^(abs(t - s))
    }
  }

  res_id <- c()
  res_t <- c()
  res_ever_treated <- c()
  res_treat <- c()
  res_X1_0 <- c()
  res_X2_0 <- c()
  res_X1 <- c()
  res_X2 <- c()
  res_y0 <- c()
  res_y <- c()

  u_mat <- Rfast::rmvnorm(N, rep(0, T), C)
  V1_mat <- Rfast::rmvnorm(N, rep(0, T), C)
  V2_mat <- Rfast::rmvnorm(N, rep(0, T), C)
  draw <- Rfast::rmvnorm(N, rep(0, 4), diag(4))
  treat_assignment <- (runif(N) > 0.5)

  for (i in 1:N) {
    id <- i
    t <- 1:T

    ever_treated <- treat_assignment[i]
    treat <- (t > T0) * ever_treated

    # factor loadings (G is for X, g is for y)
    G <- rbind(
      c(1 + draw[i, 1], draw[i, 2]),
      c(draw[i, 3], 1 + draw[i, 4])
    )
    g <- rnorm(2, mean = 1, sd = 1)

    if (parallel_trends == FALSE) {
      g <- g + c(0, 1) * ever_treated
    }

    # Generating X
    X <- (f %*% G)
    X1_0 <- X[, 1] + V1_mat[i, ]
    X2_0 <- X[, 2] + V2_mat[i, ]
    X1 <- X1_0
    # X2 <- X2_0
    X2 <- X2_0 + tau * ever_treated

    # Error term
    u <- u_mat[i, ]

    # Generating y(0)
    y0 <- X1_0 * b[1] + X2_0 * b[2] + (f %*% g) + u

    y <- Delta * treat + X1 * b[1] + X2 * b[2] + (f %*% g) + u

    res_id <- c(res_id, rep(id, T))
    res_t <- c(res_t, 1:T)
    res_ever_treated <- c(res_ever_treated, rep(ever_treated, T))
    res_treat <- c(res_treat, treat)
    res_X1_0 <- c(res_X1_0, X1_0)
    res_X2_0 <- c(res_X2_0, X2_0)
    res_X1 <- c(res_X1, X1)
    res_X2 <- c(res_X2, X2)
    res_y0 <- c(res_y0, y0)
    res_y <- c(res_y, y)
  }

  df <- data.table(
    id = res_id,
    t = res_t,
    ever_treated = res_ever_treated,
    treat = res_treat,
    X1_0 = res_X1_0,
    X2_0 = res_X2_0,
    X1 = res_X1,
    X2 = res_X2,
    y0 = res_y0,
    y = res_y
  )

  df[, rel_year := ifelse(ever_treated == TRUE, t - (T0 + 1), -Inf)]
}

estimate_ccedid <- function(df, T0 = 5) {
  df = copy(df)
  
  Fhat = df[
    ever_treated == FALSE, 
    .(ybar = mean(y), X1bar = mean(X1), X2bar = mean(X2)),
    by = "t"
  ][, 
    cbind(ybar, X1bar, X2bar)
  ]

  Fpre = Fhat[1:T0, ]
  Fpost = Fhat[(T0 +1):T, ]
  V = solve(crossprod(Fpre), t(Fpre))

  M_Fpre = diag(T0) - Fpre %*% V

  # Estimate CCEP for beta_hat
  temp = df[
    t <= T0, 
  ][,
    .(
      ytilde = as.numeric(M_Fpre %*% y), 
      X1tilde = as.numeric(M_Fpre %*% X1), 
      X2tilde = as.numeric(M_Fpre %*% X2)
    ), 
    by = "id"
  ]

  bhat = coef(feols(ytilde ~ 0 + X1tilde + X2tilde, temp, lean = TRUE))

  # imputed factor loadings and covariates
  treated_ids = df[ever_treated == TRUE, unique(id)]

  for (i in treated_ids) {
    sub = df[id == i, ]
    sub_pre = sub[t <= T0, ]

    ghat = sub_pre[, V %*% (y - cbind(X1, X2) %*% bhat)] 

    Xi_0_hat = Fhat %*% V %*% sub_pre[, cbind(X1, X2)]

    yi_0_hat = Xi_0_hat %*% bhat + Fhat %*% ghat

    df[id == i, let(
      X1_0_hat = Xi_0_hat[, 1],
      X2_0_hat = Xi_0_hat[, 2],
      y_0_hat  = as.numeric(yi_0_hat)
    )]
  }

  df[, ytilde := y - y_0_hat]
  
  # Estimate TE
  feols(ytilde ~ 0 + i(rel_year), df[ever_treated == TRUE], lean = TRUE) |> coef()
}

dfs <- lapply(1:S, \(s) {
  simulate_data(params, dgp, parallel_trends)
})
df = dfs[[1]]

# Estimate did2s, did2s w/ cov, and ccedid
ests = lapply(dfs, \(df) {
  est_did2s <- did2s:::did2s_estimate(
    df, 
    yname = "y", 
    first_stage = ~ 0 | id + t, 
    second_stage = ~ 0 + i(rel_year, ref = -Inf), 
    treatment = "treat"
  )$second_stage |> 
    coef()

  est_did2s_cov <- did2s:::did2s_estimate(
    df, 
    yname = "y", 
    first_stage = ~ X1 + X2 | id + t, 
    second_stage = ~ 0 + i(rel_year, ref = -Inf), 
    treatment = "treat"
  )$second_stage |> 
    coef()

  est_ccedid <- estimate_ccedid(df)

  return(list(
    est_did2s = est_did2s,
    est_did2s_cov = est_did2s_cov,
    est_ccedid = est_ccedid
  ))
})

T0 = params$T0
T = params$T
true_te = 
  params$tau * (1:T > T0) * (dgp >= 2) + 
  params$Delta * params$b[2] * (1:T > T0) * (dgp == 3)

bias_did2s = ests |> 
  purrr::map_dfr(\(est) est$est_did2s - true_te)
bias_did2s_cov = ests |>
  purrr::map_dfr(\(est) est$est_did2s_cov - true_te)
bias_ccedid = ests |> 
  purrr::map_dfr(\(est) est$est_ccedid - true_te)

mean_bias_did2s = colMeans(bias_did2s)
mean_bias_did2s_cov = colMeans(bias_did2s_cov)
mean_bias_ccedid = colMeans(bias_ccedid)

mse_bias_did2s = colSums(bias_did2s^2) / S
mse_bias_did2s_cov = colSums(bias_did2s_cov^2) / S
mse_bias_ccedid = colSums(bias_ccedid^2) / S

res = data.table(
  Estimator = c("TWFE", "TWFE with Covariates", "CCEDID"),
  "Mean Bias 6" = c(mean_bias_did2s[6], mean_bias_did2s_cov[6], mean_bias_ccedid[6]), 
  "MSE Bias 6" = c(mse_bias_did2s[6], mse_bias_did2s_cov[6], mse_bias_ccedid[6]), 
  "Mean Bias 7" = c(mean_bias_did2s[7], mean_bias_did2s_cov[7], mean_bias_ccedid[7]), 
  "MSE Bias 7" = c(mse_bias_did2s[7], mse_bias_did2s_cov[7], mse_bias_ccedid[7]), 
  "Mean Bias 8" = c(mean_bias_did2s[8], mean_bias_did2s_cov[8], mean_bias_ccedid[8]), 
  "MSE Bias 8" = c(mse_bias_did2s[8], mse_bias_did2s_cov[8], mse_bias_ccedid[8])
)

cols = c("Mean Bias 6", "MSE Bias 6", "Mean Bias 7", "MSE Bias 7", "Mean Bias 8", "MSE Bias 8")
res[, 
  c(cols) := lapply(.SD, \(x) round(x,2)), 
  .SDcols = cols
]

matrix_to_tex <- function(mat) {
apply(mat, 1, \(x) {
  paste0(paste0(x, collapse = " & "), " \\\\")
}) |> 
  paste0(collapse = "\n")
}

matrix_to_tex(as.matrix(res)) |>
  cat()


# gt::gt(res) |> gt::as_rtf() |> clipr::write_clip()

