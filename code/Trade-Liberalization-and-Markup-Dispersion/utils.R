#' Create matrix from vectors for rows, columns, and values
#'
#' Note the rows and columns can be anything (e.g. strings containing fips).
#' Internally, these are efficiently created to an index using
#' `indexthis::to_index` (with values 1, ..., n_unique).
#'
#' @param i Vector used for the row indices.
#' @param j Vector used for the column indices.
#' @param x Vector used for the values.
#' @param names Logical. If column and row names should be used. These
#'  correspond to the original values in `i` and `j`.
#' @param sparse Logical. Should a sparse matrix be created?
#' @example
#' fips = c("05001", "05001", "05003", "05003", "05007", "05007")
#' year = c(2010, 2011, 2010, 2011, 2010, 2011)
#' log_wage = c(2.73, 2.74, 2.8, 2.81, 2.78, 2.82)
#'
#' make_mat(fips, year, log_wage)
#' make_mat(fips, year, log_wage, names = TRUE)
#' make_mat(fips[1:5], year[1:5], log_wage[1:5], names = TRUE, sparse = TRUE)
#'
#' @return Either a matrix or a `dgCMatrix` if `sparse = TRUE`.
make_mat <- function(i, j, x = 1, names = FALSE, sparse = FALSE) {
  if (names == TRUE) {
    i_items <- indexthis::to_index(i, items = TRUE)
    i_idx <- i_items$index
    j_items <- indexthis::to_index(j, items = TRUE)
    j_idx <- j_items$index

    if (sparse == TRUE) {
      res <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x)
    } else {
      res <- matrix(NA, nrow = length(i_items$items), ncol = length(j_items$items))
      res[cbind(i_idx, j_idx)] <- x
    }
    rownames(res) <- i_items$items
    colnames(res) <- j_items$items
  } else {
    i_idx <- indexthis::to_index(i)
    j_idx <- indexthis::to_index(j)

    if (sparse == TRUE) {
      res <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x)
    } else {
      res <- matrix(NA, nrow = max(i_idx), ncol = max(j_idx))
      res[cbind(i_idx, j_idx)] <- x
    }
  }

  return(res)
}

# Copied from `punitroots::Pesaran`

#' @param TxN matrix of residuals from N single-equation models over N
#' time series of length T
#'
#' @references
#' Pesaran, M. Hashem. "General diagnostic tests for cross-sectional dependence in panels." Empirical economics 60.1 (2021): 13-50.
#'
#' @return a list containing the test-statistic (`statistic`) and the p-value (`p.value`).
#'
pesaran_test <- function(resids) {
  T <- dim(resids)[1]
  N <- dim(resids)[2]
  residCorr <- cor(resids, use = "pairwise.complete.obs")

  # test statistic
  CD <- sqrt(2 * T / N / (N - 1)) * sum(residCorr - diag(N)) / 2

  # p.value
  p.value <- 2 * (1 - pnorm(abs(CD)))

  # Output
  output <- list(
    statistic = c("CD" = CD),
    p.value = p.value
  )
  return(output)
}

#' Estimate m following Ahn and Hornstein (2013)
#' @param Z a list consisting of T \times (K+1) matrix
#'
#' @references
#' Ahn, Seung C., and Alex R. Horenstein. "Eigenvalue ratio test for the number of factors." Econometrica 81.3 (2013): 1203-1227.
#' De Vos, Ignace, Gerdie Everaert, and Vasilis Sarafidis. "A method to evaluate the rank condition for CCE estimators." Econometric Reviews 43.2-4 (2024): 123-155.
#'
DVES_m <- function(Z, m_max = ncol(Z[[1]])) {
  N <- length(Z)
  T <- nrow(Z[[1]])
  ZZt <- 1 / (N * T) * tcrossprod(do.call("cbind", Z))
  eigval <- eigen(ZZt, only.values = TRUE)$values

  h <- min(T, N * ncol(Z[[1]]))

  # could edit this to use `m_max`
  # starts at 0
  V <- function(j) { sum(eigval[(j + 1):h]) }
  GR <- sapply(1:m_max, function(j) {
    log(V(j-1) / V(j)) / log(V(j) / V(j+1))
  })
  mhat <- which.max(GR)

  return(mhat)
}


last_k_cols <- function(mat, k) {
  mat[, (ncol(mat) - k + 1):ncol(mat), drop = FALSE]
}
last_k_rows <- function(mat, k) {
  mat[(nrow(mat) - k + 1):nrow(mat), , drop = FALSE]
}

#' Estimate rho following De Vos et. al.
#' @param Z a list consisting of T \times (K+1) matrix
#' @param N the number of units averaged over
#'
#' @references
#' Robin, Jean-Marc, and Richard J. Smith. "Tests of rank." Econometric Theory 16.2 (2000): 151-175.
#' De Vos, Ignace, Gerdie Everaert, and Vasilis Sarafidis. "A method to evaluate the rank condition for CCE estimators." Econometric Reviews 43.2-4 (2024): 123-155.
#'
DVES_rho <- function(Z) {
  N <- length(Z)
  Zbar <- 1 / N * Reduce("+", Z)
  T <- nrow(Zbar)
  n <- ncol(Zbar) # K + 1
  Psi <- sqrt(T) * matrix(rnorm(T * n), nrow = n, ncol = T)
  

  # eq(13) of De Vos et al. (2024)
  # 1/N \sum_i vec(Z_i - \bar{Z}) vec(Z_i - \bar{Z})'
  Omega_i <- lapply(Z, function(Zi) {
    tcrossprod(c(last_k_rows(Zi - Zbar, n)))
  })
  # Omega_i <- lapply(Z, function(Zi) { tcrossprod(c((Zi - Zbar))) })
  Omega <- 1 / N * Reduce("+", Omega_i)

  # A = Psi Zbar Zbar' Psi'
  A <- tcrossprod(last_k_rows(Zbar, n))
  # A <- tcrossprod(Zbar)
  # B = Zbar' Psi' Psi Zbar
  B <- crossprod(last_k_rows(Zbar, n))
  # B <- crossprod(Zbar)
  eig_A <- eigen(A, symmetric = FALSE)
  eig_B <- eigen(B, symmetric = FALSE)
  # DVES_rho_seq_test(Zbar = Zbar, Omega = Omega, eig_A = eig_A, eig_B = eig_B, rho_star = 2, N = N, alpha = 0.05 / N)

  rho_star <- 0
  while (rho_star < ncol(Zbar)) {
    rho_test <- DVES_rho_seq_test(Zbar = Zbar, Omega = Omega, eig_A = eig_A, eig_B = eig_B, rho_star = rho_star, N = N, alpha = 0.05)
    # print(rho_star)
    # print(rho_test)

    if (rho_test$test_statistic > rho_test$critical_value) {
      rho_star <- rho_star + 1
    } else {
      break
    }
  }

  return(rho_star)
}


DVES_rho_seq_test <- function(Zbar, Omega, eig_A, eig_B, rho_star, N, alpha = 0.05) {
  T <- nrow(Zbar)
  n <- ncol(Zbar)
  # test-statistic in (12) of De Vos et al. (2024)
  t <- N * sum(eig_A$values[(rho_star + 1):n])

  # Weights on sum of iid \chi_1^2  r.v.
  n_eigval <- n - rho_star
  D_rho_star <- last_k_cols(eig_B$vectors, n_eigval)
  R_rho_star <- last_k_cols(eig_A$vectors, n_eigval)
  kron <- kronecker(D_rho_star, R_rho_star)
  mat <- t(kron) %*% Omega %*% kron
  w <- eigen(mat, only.values = TRUE)$values[1:n_eigval^2]

  t_bs <- sapply(1:5000, function(b) {
    sum(w * rchisq(length(w), df = 1))
  })
  return(
    list(test_statistic = t, critical_value = quantile(t_bs, 1 - alpha), alpha = alpha)
  )
}



#' Estimate rho
#' @param Z a list consisting of T \times (K+1) matrix
#' @param N the number of units averaged over
#'
#' @references
#' Robin, Jean-Marc, and Richard J. Smith. "Tests of rank." Econometric Theory 16.2 (2000): 151-175.
#' De Vos, Ignace, Gerdie Everaert, and Vasilis Sarafidis. "A method to evaluate the rank condition for CCE estimators." Econometric Reviews 43.2-4 (2024): 123-155.
#'
rho_test <- function(Z) {
  N <- length(Z)
  kappa <- N^(-1 / 4)

  # CCE cross-sectional averages \hat{P}_n \to^p P_0 = CF
  Zbar <- 1 / N * Reduce("+", Z)
  m <- nrow(Zbar)
  k <- ncol(Zbar)

  # eq(13) of De Vos et al. (2024)
  # 1/N \sum_i vec(Z_i - \bar{Z}) vec(Z_i - \bar{Z})'
  # vec(cal(M)) ~ N(0, Omega)
  Omega <- 1 / N * Reduce("+", lapply(Z, \(Zi) tcrossprod(c(Zi - Zbar))))

  # SVD-Decomposition of Zbar
  svd_Zbar <- svd(Zbar)
  Sigma <- Matrix::Diagonal(length(svd_Zbar$d), svd_Zbar$d)
  P <- svd_Zbar$u
  Q <- svd_Zbar$v
  # P %*% Sigma %*% t(Q)

  r_hat <- sum(svd_Zbar$d >= kappa)

  P_0_2 <- ncol(P) - r_0

  # \Pi_n = \sum_i Z_i where Z_i is a T x K+1 matrix (m x k in their notation)
  # test-statistic first sentence of pg 1794 of Chen and Fang (2019) with \tau_N = sqrt(N)
  idx_tail <- (rho_star + 1):length(eigval)
  t <- N * sum(eigval[idx_tail]^2)

  tb <- sapply(1:100, function(b) {
    Mb <- matrix(mvtnorm::rmvnorm(1, sigma = Omega), nrow = nrow(Zbar), ncol = ncol(Zbar))
  })
  # for (b in 1:B) {
  #
  # }





  # n - rho_sta smallest eigen-vectors
  D_rho_star <- eig_tZZ$vectors[, (rho_star + 1):n, drop = FALSE]
  R_rho_star <- eig_ZtZ$vectors[, (rho_star + 1):n, drop = FALSE]
  kron <- kronecker(D_rho_star, R_rho_star)
  crossprod(kron, Omega) %*% kron

  # Simulate weighted-sum of independent \chi^2_1
}
