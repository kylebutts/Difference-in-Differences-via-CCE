#' Draw from the AR(1) distribution
#' @param T Integer number of observations to draw
#' @param rho Vector of autocorrelation coefficients
rAR <- function(T, rho = 0.75) {
  x <- arima.sim(n = T, list(order = c(1, 0, 0), ar = rho), n.start = 10)
  return(as.numeric(x))
}

# %%
convert_to_data_frame <- function(y, X, d_it) {
  N <- ncol(y)
  T <- nrow(y)
  df <- data.frame(
    id = rep(1:N, each = T),
    t = rep(1:T, times = N),
    y = c(y),
    X = c(X),
    treat = c(d_it)
  )
}


# %%
library(fixest)
est_twfe_imputation <- function(y, X, d_it) {
  df <- convert_to_data_frame(y, X, d_it)
  fs <- feols(y ~ 1 | id + t, df[df$treat == 0, ])
  df$y_resid = df$y - predict(fs, newdata = df)
  ss <- feols(y_resid ~ i(treat), data = df)
  est <- unname(coef(ss)["treat::1"])

  return(est)
}
est_twfe_imputation_cov <- function(y, X, d_it) {
  df <- convert_to_data_frame(y, X, d_it)
  fs <- feols(y ~ X | id + t, df[df$treat == 0, ])
  df$y_resid = df$y - predict(fs, newdata = df)
  ss <- feols(y_resid ~ i(treat), data = df)
  est <- unname(coef(ss)["treat::1"])

  return(est)
}

# Confirm
# df <- convert_to_data_frame(y, X, d_it)
# solve(crossprod(cbind(1, c(d_it), c(X))), crossprod(cbind(1, c(d_it), c(X)), c(y)))
# fixest::feols(y ~ 1 + i(treat) + X, df)

# %%
est_tecce <- function(y, X, d_it) {
  D <- colSums(d_it) > 0
  T_0 <- max(which(rowSums(d_it) == 0))
  T <- nrow(y)
  N <- ncol(y)

  Fhat <- cbind(
    rowMeans(y[, which(!D), drop = FALSE]),
    rowMeans(X[, which(!D), drop = FALSE])
  )
  Fhat_pre <- Fhat[1:T_0, , drop = FALSE]

  y0_hat <- matrix(NA, nrow = T, ncol = N)
  ai_solver <- solve(crossprod(Fhat_pre), t(Fhat_pre))
  for (i in which(D)) {
    ai_hat <- ai_solver %*% y[1:T_0, i]
    y0_hat[, i] <- Fhat %*% ai_hat
  }

  di <- y - y0_hat
  est <- mean(di[(T_0 + 1):T, which(D)])

  return(est)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Xu by hand
est_xu_with_X <- function(y, X, d_it, r = 1) {
  T <- nrow(y)
  N <- ncol(y)
  k <- ncol(X) / N

  # Bai (2009) procedure to estimate factors and loadings using control group
  control_idx <- which(colSums(d_it) == 0)
  y_control <- y[, control_idx]
  X_control <- X[, control_idx]

  ## Initial estimate of beta, factors, and loadings
  vy <- as.vector(y_control)
  ## Temp: subtract grand mean
  ## vy <- vy - mean(vy)
  vx <- as.vector(X_control)
  ## Temp: subtract grand mean
  ## vx <- vx - mean(vx)
  inv_xx <- solve(crossprod(vx, vx))
  b_0 <- inv_xx %*% crossprod(vx, vy)
  w <- matrix(vy - vx %*% b_0, ncol = length(control_idx), byrow = FALSE)

  svd_res <- svd(w %*% t(w))
  f <- svd_res$u[, 1:r, drop = FALSE] * sqrt(T)
  l <- t(w) %*% f / T
  b <- inv_xx %*% crossprod(vx, vy - as.vector(f %*% t(l)))

  ## Iterate until convergence
  j <- 1
  while (max(abs(b - b_0)) > 1e-3 && j <= 1000) {
    b_0 <- b
    w <- matrix(vy - vx %*% b_0, ncol = length(control_idx), byrow = FALSE)
    svd_res <- svd(w %*% t(w))
    f <- svd_res$u[, 1:r, drop = FALSE] * sqrt(T)
    l <- t(w) %*% f / T
    b <- inv_xx %*% crossprod(vx, vy - as.vector(f %*% t(l)))
    j <- j + 1
  }

  # Estimate treatment effects
  T_0 <- max(which(rowSums(d_it) == 0))
  f_0 <- f[1:T_0, , drop = FALSE]
  treated_idx <- which(colSums(d_it) > 0)

  yh <- matrix(NA, nrow = T, ncol = length(treated_idx))
  for (i in seq_along(treated_idx)) {
    unit_idx <- treated_idx[i]
    X_unit <- X[, (unit_idx - 1) * k + (1:k), drop = FALSE]
    ai <- solve(
      crossprod(f_0),
      crossprod(
        f_0,
        y[1:T_0, unit_idx] -
          X[1:T_0, (unit_idx - 1) * k + (1:k), drop = FALSE] %*% b
      )
    )
    yh[, i] <- X_unit %*% b + f %*% ai
  }

  di <- y[, treated_idx] - yh
  d <- rowMeans(di)

  # Overall ATT (average of post-treatment effects)
  post_period <- (T_0 + 1):T
  overall_att <- mean(d[post_period])

  return(overall_att)
}

est_xu <- function(y, d_it, r = 1) {
  T <- nrow(y)
  N <- ncol(y)

  # Bai (2009) procedure to estimate factors and loadings using control group
  control_idx <- which(colSums(d_it) == 0)
  y_control <- y[, control_idx]

  ## Initial estimate of factors and loadings
  vy <- as.vector(y_control)
  w <- matrix(vy, ncol = length(control_idx), byrow = FALSE)

  svd_res <- svd(w %*% t(w))
  f <- svd_res$u[, 1:r, drop = FALSE] * sqrt(T)
  l <- t(w) %*% f / T

  # Estimate treatment effects
  T_0 <- max(which(rowSums(d_it) == 0))
  f_0 <- f[1:T_0, , drop = FALSE]
  treated_idx <- which(colSums(d_it) > 0)

  yh <- matrix(NA, nrow = T, ncol = length(treated_idx))
  for (i in seq_along(treated_idx)) {
    unit_idx <- treated_idx[i]
    ai <- solve(
      crossprod(f_0),
      crossprod(
        f_0,
        y[1:T_0, unit_idx]
      )
    )
    yh[, i] <- f %*% ai
  }

  di <- y[, treated_idx] - yh
  d <- rowMeans(di)

  # Overall ATT (average of post-treatment effects)
  post_period <- (T_0 + 1):T
  overall_att <- mean(d[post_period])

  return(overall_att)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_xu <- function(y, X, d_it, r = 1) {
  T <- nrow(y)
  N <- ncol(y)
  k <- ncol(X) / N

  # Bai (2009) procedure to estimate factors and loadings using control group
  control_idx <- which(colSums(d_it) == 0)
  y_control <- y[, control_idx]
  X_control <- X[, control_idx]

  ## Initial estimate of beta, factors, and loadings
  vy <- as.vector(y_control)
  ## Temp: subtract grand mean
  ## vy <- vy - mean(vy)
  vx <- as.vector(X_control)
  ## Temp: subtract grand mean
  ## vx <- vx - mean(vx)

  inv_xx <- solve(crossprod(vx, vx))
  b_0 <- inv_xx %*% crossprod(vx, vy)
  w <- matrix(vy - vx %*% b_0, ncol = length(control_idx), byrow = FALSE)

  svd_res <- svd(w %*% t(w))
  f <- svd_res$u[, 1:r, drop = FALSE] * sqrt(T)
  l <- t(w) %*% f / T
  b <- inv_xx %*% crossprod(vx, vy - as.vector(f %*% t(l)))

  ## Iterate until convergence
  j <- 1
  while (max(abs(b - b_0)) > 1e-3 && j <= 1000) {
    b_0 <- b
    w <- matrix(vy - vx %*% b_0, ncol = length(control_idx), byrow = FALSE)
    svd_res <- svd(w %*% t(w))
    f <- svd_res$u[, 1:r, drop = FALSE] * sqrt(T)
    l <- t(w) %*% f / T
    b <- inv_xx %*% crossprod(vx, vy - as.vector(f %*% t(l)))
    j <- j + 1
  }

  # Estimate treatment effects
  T_0 <- max(which(rowSums(d_it) == 0))
  f_0 <- f[1:T_0, , drop = FALSE]
  treated_idx <- which(colSums(d_it) > 0)

  yh <- matrix(NA, nrow = T, ncol = length(treated_idx))
  for (i in seq_along(treated_idx)) {
    unit_idx <- treated_idx[i]
    X_unit <- X[, (unit_idx - 1) * k + (1:k), drop = FALSE]
    ai <- solve(
      crossprod(f_0),
      crossprod(
        f_0,
        y[1:T_0, unit_idx] -
          X[1:T_0, (unit_idx - 1) * k + (1:k), drop = FALSE] %*% b
      )
    )
    yh[, i] <- X_unit %*% b + f %*% ai
  }

  di <- y[, treated_idx] - yh
  d <- rowMeans(di)

  # Overall ATT (average of post-treatment effects)
  post_period <- (T_0 + 1):T
  overall_att <- mean(d[post_period])

  Y.ct.bar <- rowMeans(yh)

  return(list(beta_xu = b, Y.ct.bar_xu = Y.ct.bar))
}

test_gsynth_with_x <- function(y, X, d_it, r = 1) {
  df <- convert_to_data_frame(y, X, d_it)
  res <- gsynth::gsynth(
    data = df,
    Y = "y",
    D = "treat",
    X = "X",
    index = c("id", "t"),
    force = "none",
    estimator = "ife",
    r = r,
    CV = FALSE,
    min.T0 = 4,
    se = FALSE
  )

  return(list(
    beta_gsynth = res$beta,
    Y.ct.bar_gsynth = res$Y.bar[, "Y.ct.bar"]
  ))
}


# %%
#' @param n Number of units
#' @param t Number of time periods
#' @param fact_label String containing the factor label specification. One of:
#'   - "Constant"
#'   - "AR(1) Process"
#'   - "Linear Trend"
#' @param loading_shift Either 0, -1, or 1. Level-shift in factor-loading for treated group. Induces non-parallel trend when fact = 2 or 3
#' @param d_0x Treatment effect on x
#' @param d_0y Treatment effect on y
#' @param verbose Print extra information
run_simulations <- function(
  nsim = 50,
  N = 100,
  T = 8,
  T_0 = T - 1,
  dgp_number = 1,
  fact_label = c("Constant", "Linear Trend", "Autocorrelated"),
  eta_y = 1,
  tau_x = 0,
  avg_diff_in_loading = 0,
  eps_autocorrelated = FALSE
) {
  #
  beta <- 1
  true_att <- eta_y + beta * tau_x

  # Generate treatment dummy (half get treated)
  D_i <- c(rep(0, N / 2), rep(1, N / 2))
  post <- as.numeric((1:T) > T_0)

  # T x N matrix = D_i * post_t
  d_it <- outer(post, D_i)

  # First, generate a single factor
  # (will be held constant across simulations)
  F <- NULL
  if (fact_label == "Constant") {
    F <- rep(1, 15)
  } else if (fact_label == "Linear Trend") {
    F <- 1 + (1:15 / 8)
  } else if (fact_label == "Autocorrelated") {
    F <- 1.5 + rAR(15, rho = 0.6)
  }
  ## Keep the last 10 elements so that F_T are the same across values of T
  F <- cbind(F[(15 - T + 1):15])

  bias1 <- rep(0, nsim)
  bias2 <- rep(0, nsim)
  bias3 <- rep(0, nsim)
  bias4 <- rep(0, nsim)
  bias5 <- rep(0, nsim)
  bias6 <- rep(0, nsim)
  for (is in 1:nsim) {
    # Generate factor-loadings
    # Note \alpha_i and \lambda_i are drawn to be correlated with one another
    alpha_i <- vector("numeric", length = N)
    lambda_i <- vector("numeric", length = N)
    for (i in 1:N) {
      mu <- 2 + (D_i[i] * avg_diff_in_loading)
      draw <- MASS::mvrnorm(
        n = 1,
        mu = rep(mu, 2),
        Sigma = cbind(c(0.5, 0.25), c(0.25, 0.5))
      )
      alpha_i[i] <- draw[1]
      lambda_i[i] <- draw[2]
    }

    # Draw disturbance terms
    if (eps_autocorrelated == TRUE) {
      eps_y <- unlist(lapply(1:N, function(i) {
        rAR(T, rho = 0.4)
      }))
      eps_x <- unlist(lapply(1:N, function(i) {
        rAR(T, rho = 0.4)
      }))
    } else {
      eps_y <- matrix(rnorm(T * N, mean = 0, sd = 0.4), nrow = T, ncol = N)
      eps_x <- matrix(rnorm(T * N, mean = 0, sd = 0.4), nrow = T, ncol = N)
    }

    # T X N matrix of X_it
    X <- F %*% t(lambda_i) + d_it * tau_x + eps_x

    # T X N matrix of y_it
    y <- X * beta + F %*% t(alpha_i) + d_it * eta_y + eps_y
    y0 <- X * beta + F %*% t(alpha_i) + eps_y

    tau_hat_1 <- est_tecce(y, X, d_it)
    tau_hat_2 <- est_xu(y, d_it, r = 1)
    tau_hat_3 <- est_xu_with_X(y, X, d_it, r = 1)
    tau_hat_4 <- est_twfe_imputation(y, X, d_it)
    tau_hat_5 <- est_twfe_imputation_cov(y, X, d_it)

    bias1[is] <- (tau_hat_1 - true_att)
    bias2[is] <- (tau_hat_2 - true_att)
    bias3[is] <- (tau_hat_3 - true_att)
    bias4[is] <- (tau_hat_4 - true_att)
    bias5[is] <- (tau_hat_5 - true_att)

    ## Testing
    # res_test_xu <- test_xu(y, X, d_it, r = 1)
    # cat("xu: \n")
    # print(res_test_xu)
    # cat("\n")
    #
    # res_test_gsynth <- test_gsynth_with_x(y, X, d_it, r = 1)
    # cat("gsynth: \n")
    # print(res_test_gsynth)
    # cat("\n")
    #
    # res_test_gsynth_r_2 <- test_gsynth_with_x(y, X, d_it, r = 2)
    # cat("gsynth: \n")
    # print(res_test_gsynth)
    # cat("\n")
    #
    # ggplot() +
    # geom_point(aes(x = 1:T, y = res_test_xu$Y.ct.bar_xu, color = "By hand")) +
    # geom_point(aes(
    # x = 1:T,
    # y = res_test_gsynth$Y.ct.bar_gsynth,
    # color = "gsynth"
    # )) +
    # geom_point(aes(
    # x = 1:T,
    # y = res_test_gsynth_r_2$Y.ct.bar_gsynth,
    # color = "gsynth (r = 2)"
    # )) +
    # geom_point(aes(
    # x = 1:T,
    # y = rowMeans(y0[, D_i == 1]),
    # color = "true trend for treated"
    # )) +
    # labs(x = NULL, color = NULL, y = "\\bar{Y}(0) for treated units") +
    # kfbmisc::theme_kyle(base_size = 14, legend = "top")
  }

  mean_biases <- c(
    mean(bias1),
    mean(bias2),
    mean(bias3),
    mean(bias4),
    mean(bias5)
  )
  mean_abs_biases <- c(
    mean(abs(bias1)),
    mean(abs(bias2)),
    mean(abs(bias3)),
    mean(abs(bias4)),
    mean(abs(bias5))
  )
  mean_rmses <- c(
    sqrt(mean(bias1^2)),
    sqrt(mean(bias2^2)),
    sqrt(mean(bias3^2)),
    sqrt(mean(bias4^2)),
    sqrt(mean(bias5^2))
  )
  bias_emp_ci_lower <- unname(c(
    quantile(bias1, 0.025),
    quantile(bias2, 0.025),
    quantile(bias3, 0.025),
    quantile(bias4, 0.025),
    quantile(bias5, 0.025)
  ))
  bias_emp_ci_upper <- unname(c(
    quantile(bias1, 0.975),
    quantile(bias2, 0.975),
    quantile(bias3, 0.975),
    quantile(bias4, 0.975),
    quantile(bias5, 0.975)
  ))

  ## Asymptotic confidence interval vs. empirical
  ## mean_biases - 1.96 * mean_rmses
  ## bias_emp_ci_lower

  res <- tibble(
    estimator = c(
      "TECCE",
      "gsynth",
      "gsynth with Control",
      # "xu",
      "TWFE via OLS",
      "TWFE with Control"
    ),
    bias = mean_biases,
    rmse = mean_rmses,
    bias_emp_ci_lower = bias_emp_ci_lower,
    bias_emp_ci_upper = bias_emp_ci_upper,
    abs_bias = mean_abs_biases,
  )

  ## Params from the simulation
  res$dgp_number = dgp_number
  res$N = N
  res$T = T
  res$fact_label = fact_label
  res$eps_autocorrelated = eps_autocorrelated
  res$avg_diff_in_loading = avg_diff_in_loading
  res$eta_y = eta_y
  res$tau_x = tau_x
  res$beta = beta
  res$true_att = true_att

  return(res)
}
