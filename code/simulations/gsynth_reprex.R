## This code highlights a subtle bug in `gsynth` involving the grand-demeaning when `force = "none"` and `r = 1`.
library(gsynth)

N = 100
T = 8
B = 250

atts = vector("numeric", B)
atts_r_2 = vector("numeric", B)
atts_unit_fe = vector("numeric", B)
for (b in 1:B) {
  D_i <- c(rep(0, N / 2), rep(1, N / 2))
  post <- +(1:T > (T - 1))
  d_it <- outer(post, D_i)

  # Generate factors and factor-loadings
  F <- cbind(1:T)
  alpha_i <- rnorm(N, 1 + 1 * D_i, sd = 0.5)

  # T X N matrix of y_it with no treatment effects
  eps_y <- matrix(rnorm(T * N, mean = 0, sd = 0.4), nrow = T, ncol = N)
  y <- F %*% t(alpha_i) + eps_y

  df <- data.frame(
    id = rep(1:N, each = T),
    t = rep(1:T, times = N),
    y = c(y),
    treat = c(d_it)
  )

  ## Default correct choice of r
  res <- gsynth(
    y ~ treat,
    data = df,
    index = c("id", "t"),
    force = "none", ## no unit or time FE
    r = 1,
    CV = FALSE ## 1 factor
  )
  atts[b] <- res$att.avg

  ## Adding second factor
  res <- gsynth(
    y ~ treat,
    data = df,
    index = c("id", "t"),
    force = "none", ## no unit or time FE
    r = 2,
    CV = FALSE ## 2 factors
  )
  atts_r_2[b] <- res$att.avg

  ## Adding unit FE
  res <- gsynth(
    y ~ treat,
    data = df,
    index = c("id", "t"),
    force = "unit",
    r = 1,
    CV = FALSE ## 1 factor
  )
  atts_unit_fe[b] <- res$att.avg
}

hist(atts)
# hist(atts_r_2)
# hist(atts_unit_fe)
