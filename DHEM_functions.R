## ================================================================
## DHEM_functions.R
##
## All estimation / simulation functions used to produce the figures
## and tables in the paper are collected here, copied unmodified from
## their original source files:
##   - gmm_functions.R  -> Gaussian Mixture Model  (prefix: gmm_)
##   - wmm_functions.R  -> Weibull Mixture Model    (prefix: wmm_ / weibull_ / barrier_)
##   - ZIP_functions.R  -> Zero-Inflated Poisson    (prefix: zip_)
##
## Function bodies are copied as-is; only these section banners were
## added so the file can be sourced as a single entry point from
## DHEM.R. Do not edit the functions themselves when reproducing the
## paper's results.
## ================================================================

## ----------------------------------------------------------------
## Section 1: Gaussian Mixture Model (GMM) functions
## Source: gmm_functions.R
## ----------------------------------------------------------------
library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)


## ------------------------------------------------------------
## Data generator + init
## ------------------------------------------------------------

invSigma_times <- function(R, v) {
  # v: length d
  # --- likely error point: dimension mismatch between v and R ---
  y <- backsolve(R, v, transpose = TRUE)  # solves R^T y = v
  backsolve(R, y, transpose = FALSE)      # solves R z = y  -> z = Sigma^{-1} v
}

rgmm <- function(n, pi, mu_list, Sigma_list) {
  K <- length(pi)
  d <- length(mu_list[[1]])
  z <- sample.int(K, size=n, replace=TRUE, prob=pi)
  x <- matrix(0, n, d)
  for (k in 1:K) {
    idx <- which(z == k)
    if (!length(idx)) next
    cholS <- chol(Sigma_list[[k]])
    Z <- matrix(rnorm(length(idx)*d), length(idx), d)
    x[idx, ] <- sweep(Z %*% cholS, 2, mu_list[[k]], "+")
  }
  list(x=x)
}

# Initialize GMM parameters: generate pi, mu, and Sigma with minimum separation constraint on means
gmm_init <- function(x, K,
                     pi_min = 0.05,
                     delta = 1,
                     Sigma_scale = 3,
                     max_tries = 5000L,
                     seed = NULL
){
  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  if (K < 1) stop("K must be >= 1")
  if (pi_min < 0) stop("pi_min must be >= 0")
  if (K * pi_min > 1 + 1e-12) stop("Infeasible: K * pi_min must be <= 1")
  if (!is.finite(delta) || delta <= 0) stop("delta must be > 0")
  if (!is.finite(Sigma_scale) || Sigma_scale <= 0) stop("Sigma_scale must be > 0")

  # Sample covariance matrix
  S <- tryCatch(cov(x), error = function(e) diag(d))

  ## 1) pi: original scheme (ensure minimum pi_min for each component,
  ##    distribute the remaining mass randomly)
  if (K == 1) {
    pi <- 1
  } else {
    rem <- 1 - K * pi_min
    w <- rgamma(K, shape = 1, rate = 1)
    w <- w / sum(w)
    pi <- pi_min + rem * w
  }

  ## 2) Sigma: replicate the sample covariance matrix S for all K components
  Sigma_list <- replicate(K, S, simplify = FALSE)

  ## 3) mu: generate K points from N(mean(x), Sigma_scale * S),
  ##    repeat until all pairwise distances are at least delta
  xbar <- colMeans(x)
  Sig_mu <- Sigma_scale * S

  rmvnorm_one <- function(mu, Sigma) {
    # Base R only: sampling via Cholesky decomposition
    R <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(R)) {
      # If not numerically PSD, add diagonal jitter
      eps <- 1e-8
      Sigma2 <- Sigma + diag(eps, d)
      R <- chol(Sigma2)
    }
    as.numeric(mu + drop(t(R) %*% rnorm(d)))
  }

  mu_list <- vector("list", K)
  mu_list[[1]] <- as.numeric(xbar)

  for (k in 2:K) {
    ok <- FALSE
    for (t in 1:max_tries) {
      cand <- rmvnorm_one(xbar, Sig_mu)
      # Euclidean distance criterion: at least delta away from all previous mu's
      if (all(sapply(1:(k - 1), function(j) sum((cand - mu_list[[j]])^2) >= delta/10))) {
        mu_list[[k]] <- cand
        ok <- TRUE
        break
      }
    }
    if (!ok) stop(sprintf("Failed to generate mu with pairwise distance >= delta after %d tries (k=%d).",
                          max_tries, k))
  }

  list(pi = pi, mu = mu_list, Sigma = Sigma_list)
}


make_data_and_init <- function(n, K, theta_true, delta,seed=sample.int(1e9, 1)) {
  set.seed(seed)
  # cat("--- ",theta_true)
  x <- rgmm(n, theta_true$pi, theta_true$mu, theta_true$Sigma)$x

  set.seed(seed)
  theta_init <- gmm_init(x, K, delta = delta * 1.1)

  list(x = x, theta_init = theta_init,seed=seed)
}


dmvnorm_log <- function(X, m, S) {
  # X: n x d matrix
  # m: length d vector
  # S: d x d SPD matrix
  # return: length n vector of log N(x_i | m, S)

  X <- as.matrix(X)
  m <- as.numeric(m)
  d <- ncol(X)

  R <- chol(S)
  Xm <- sweep(X, 2, m, "-")                 # n x d
  Y  <- backsolve(R, t(Xm), transpose = TRUE) # d x n, solves R^T y = (x-m)^T
  quad <- colSums(Y^2)                      # length n
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (d * log(2 * pi) + logdet + quad)
}

gmm_Qr <- function(x, pi, mu, sigma, latent) {
  x <- as.matrix(x)
  n <- nrow(x)
  K <- length(pi)

  Q <- 0
  for (k in 1:K) {
    # log N(x_i | mu_k, Sigma_k) :n vector
    log_comp <- dmvnorm_log(x, mu[[k]], sigma[[k]])
    Q <- Q + sum(latent[, k] * (log(pi[k]) + log_comp))
  }
  Q
}



sumlog_mahalanobis_barrier_k <- function(k, mu, sigma, delta) {
  # b_k(mu_k) = log( sum_{l!=k} (mu_k-mu_l)' Sigma_k^{-1} (mu_k-mu_l) - delta )
  K <- length(mu)
  mu_k <- as.numeric(mu[[k]])
  Sig_k <- sigma[[k]]

  Rk <- chol(Sig_k)
  s <- 0

  for (l in 1:K) {
    if (l == k) next
    diff <- as.numeric(mu_k - mu[[l]])
    yk <- backsolve(Rk, diff, transpose = TRUE)  # Rk^T y = diff
    d2 <- sum(yk^2)                              # diff' Sig_k^{-1} diff
    s <- s + d2
  }

  inside <- s - delta
  if (inside <= 0) return(-Inf)
  log(inside)
}


gmm_BQr<- function(x, pi, mu, sigma, latent, bw, delta = 0.1) {
  x <- as.matrix(x)
  K <- length(pi)

  BQ <- 0
  for (k in 1:K) {
    log_comp <- dmvnorm_log(x, mu[[k]], sigma[[k]])
    BQ <- BQ + sum(latent[, k] * (log(pi[k]) + log_comp))

    bk <- sumlog_mahalanobis_barrier_k(k = k, mu = mu, sigma = sigma, delta = delta)
    if (!is.finite(bk)) return(-Inf)
    BQ <- BQ + bw * bk
  }
  BQ
}


solve_BQ_sumlog_mubarrier_newton <- function(k, x, mu, sigma, latent,
                                             delta, bw,
                                             maxiter = 100, tol = 1e-6,
                                             alpha_shrink = 0.5
) {
  x <- as.matrix(x)
  d <- ncol(x)
  K <- length(mu)

  # weights
  w <- latent[, k]
  Nk <- sum(w)
  if (Nk <= 0) return(as.numeric(mu[[k]]))

  s <- colSums(x * w)  # s_k

  # fixed Sigma_k and its chol
  Sig_k <- sigma[[k]]
  Rk <- chol(Sig_k)

  A_times <- function(v) invSigma_times(Rk, v)  # Sigma_k^{-1} v

  # other means fixed at this inner solve
  mu_others <- mu[-k]
  m <- Reduce(`+`, lapply(mu_others, as.numeric))
  Km1 <- K - 1L

  u_fun <- function(mk) Km1 * mk - m  # sum_{l!=k} (mk - mu_l)

  D_fun <- function(mk) {
    # D(mk) = sum_{l!=k} (mk-mu_l)' Sig_k^{-1} (mk-mu_l) - delta
    acc <- 0
    for (ml in mu_others) {
      diff <- mk - as.numeric(ml)
      Ad <- A_times(diff)
      acc <- acc + sum(diff * Ad)
    }
    acc - delta
  }

  G_fun <- function(mk) {
    u <- u_fun(mk)
    D <- D_fun(mk)
    D * (s - Nk * mk) + 2 * bw * u
  }

  J_fun <- function(mk) {
    # J_G(mk) = (-Nk*D + 2*bw*(K-1)) I + 2 (s - Nk*mk) (A u)^T
    u <- u_fun(mk)
    D <- D_fun(mk)
    Au <- A_times(u)

    base <- (-Nk * D + 2 * bw * Km1)
    base * diag(d) + 2 * ((s - Nk * mk) %*% t(Au))
  }

  mk <- as.numeric(mu[[k]])

  # Ensure initial feasibility: if not feasible, do small damping toward weighted mean
  if (D_fun(mk) <= 0) {
    mk_em <- s / Nk
    # try convex combination toward mk_em until feasible
    alpha <- 1
    for (bt in 1:30) {
      mk_try <- (1 - alpha) * mk + alpha * mk_em
      if (D_fun(mk_try) > 0) { mk <- mk_try; break }
      alpha <- alpha * 0.5
    }
    if (D_fun(mk) <= 0) return(as.numeric(mu[[k]]))  # fail-safe: no move
  }

  # Damped Newton iterations
  for (it in 1:maxiter) {
    G <- G_fun(mk)
    if (!all(is.finite(G))) break
    if (sqrt(sum(G^2)) < tol) break

    J <- J_fun(mk)
    if (!all(is.finite(J))) break

    step <- tryCatch(solve(J, -G), error = function(e) rep(0, d))

    # backtracking: keep D>0 and reduce ||G||
    alpha <- 1
    normG <- sqrt(sum(G^2))
    repeat {
      mk_try <- mk + alpha * step
      D_try <- D_fun(mk_try)
      if (D_try > 0) {
        G_try <- G_fun(mk_try)
        if (all(is.finite(G_try)) && sqrt(sum(G_try^2)) < normG) break
      }
      alpha <- alpha * alpha_shrink
      if (alpha < 1e-10) { mk_try <- mk; break }
    }

    mk <- mk_try
  }

  mk
}


is_feasible_sym <- function(mu, sigma, delta) {
  K <- length(mu)

  for (k in 1:(K - 1)) for (l in (k + 1):K) {
    diff <- as.numeric(mu[[k]] - mu[[l]])

    # --- likely error point: sigma[[k]] not SPD -> chol fails ---
    Rk <- chol(sigma[[k]])
    yk <- backsolve(Rk, diff, transpose = TRUE)
    d2k <- sum(yk^2)

    # --- likely error point: sigma[[l]] not SPD -> chol fails ---
    Rl <- chol(sigma[[l]])
    yl <- backsolve(Rl, diff, transpose = TRUE)
    d2l <- sum(yl^2)

    if ((d2k + d2l) <= delta) return(FALSE)
  }
  TRUE
}

is_feasible_sumlog <- function(mu, sigma, delta) {
  K <- length(mu)
  for (k in 1:K) {
    val <- sumlog_mahalanobis_barrier_k(k, mu, sigma, delta)
    if (!is.finite(val)) return(FALSE)
  }
  TRUE
}

# Accept a proposed (mu, Sigma) update only if it is feasible and does not decrease the barrier-augmented objective
accept_mu_sigma_by_damping_feasible <- function(mu_old, sigma_old,
                                                x, pi_new, latent,
                                                mu_prop,
                                                bw, delta,
                                                BQ_old,
                                                alpha0 = 1.0,
                                                alpha_shrink = 0.5,
                                                max_backtrack = 20
) {
  K <- length(mu_old)

  # compute Sigma and objective for the proposal
  sigma_prop <- gmm_sigma_update(x, mu_prop, latent)

  # check feasibility of the proposal
  feasible_prop <- is_feasible_sym(mu_prop, sigma_prop, delta)

  if (feasible_prop) {
    BQ_prop <- gmm_BQr(x, pi_new, mu_prop, sigma_prop,
                       latent = latent, bw = bw, delta = delta)

    # accept immediately if feasible and non-decreasing
    if (is.finite(BQ_old) && is.finite(BQ_prop) && (BQ_prop >= BQ_old)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }

    # accept if the old objective is not finite but the proposal is finite
    if (!is.finite(BQ_old) && is.finite(BQ_prop)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
  }

  # backtracking: interpolate between old and proposed mu
  alpha <- alpha0
  mu_try <- mu_old
  sigma_try <- sigma_old
  BQ_try <- BQ_old

  for (bt in 1:max_backtrack) {
    alpha <- alpha * alpha_shrink

    mu_try <- lapply(seq_len(K), function(k) {
      (1 - alpha) * as.numeric(mu_old[[k]]) + alpha * as.numeric(mu_prop[[k]])
    })

    sigma_try <- gmm_sigma_update(x, mu_try, latent)

    # skip infeasible candidates
    if (!is_feasible_sym(mu_try, sigma_try, delta)) next

    BQ_try <- gmm_BQr(x, pi_new, mu_try, sigma_try,
                      latent = latent, bw = bw, delta = delta)

    if (is.finite(BQ_try) && is.finite(BQ_old) && (BQ_try >= BQ_old)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }

    if (!is.finite(BQ_old) && is.finite(BQ_try)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
  }

  # fallback: keep the old parameters
  return(list(mu = mu_old, sigma = sigma_old, BQ = BQ_old, alpha = 0.0))
}


gmm_B_value <- function(theta, delta = 0.1, symmetric = TRUE) {
  mu    <- theta$mu
  sigma <- theta$Sigma
  K <- length(mu)

  B <- 0
  for (k in seq_len(K)) {
    bk = sumlog_mahalanobis_barrier_k(k = k, mu = mu, sigma = sigma, delta = delta)
    B <- B + bk
  }
  B
}

# Delta |B(theta0, theta1)| = | B(theta1) - B(theta0) |
gmm_delta_abs_B <- function(theta0, theta1, delta = 0.1, symmetric = TRUE) {
  B0 <- gmm_B_value(theta0, delta = delta, symmetric = symmetric)
  B1 <- gmm_B_value(theta1, delta = delta, symmetric = symmetric)

  if (!is.finite(B0) || !is.finite(B1)) return(Inf)

  abs(B1 - B0)
}


# Cap the barrier weight so that the barrier change does not exceed the KL-based target
cap_bw_by_barrier <- function(theta0, theta1, eta, dkl_1, bw, delta, symmetric = TRUE) {
  absDeltaB <- gmm_delta_abs_B(theta0, theta1, delta = delta, symmetric = symmetric)

  # no restriction if the barrier change is zero or not finite
  if (!is.finite(absDeltaB) || absDeltaB <= 0) {
    return(list(bw_new = bw, violated = FALSE, absDeltaB = absDeltaB))
  }

  # delta_target <- eta * dkl_1
  delta_target <- 0.5*eta * dkl_1

  bw_cap <- delta_target / absDeltaB
  bw_new <- min(bw, bw_cap)

  violated <- (bw > bw_cap)
  list(bw_new = bw_new, violated = violated, absDeltaB = absDeltaB, bw_cap = bw_cap, delta_target = delta_target)
}
# ------------------------------------------------------------------
# Sumlog-barrier version of accept_mu_sigma_by_damping_feasible()
#   - feasibility: is_feasible_sumlog()
#   - objective  : gmm_BQr()
# ------------------------------------------------------------------
accept_mu_sigma_by_damping_feasible_sumlog <- function(mu_old, sigma_old,
                                                       x, pi_new, latent,
                                                       mu_prop,
                                                       bw, delta,
                                                       BQ_old,
                                                       alpha0 = 1.0,
                                                       alpha_shrink = 0.5,
                                                       max_backtrack = 20) {
  K <- length(mu_old)

  # 1) proposal sigma + feasibility + BQ
  sigma_prop <- gmm_sigma_update(x, mu_prop, latent)

  # sumlog feasibility (all k: sum_{l!=k} d2_{kl}(Sigma_k^{-1}) - delta > 0)
  feasible_prop <- is_feasible_sumlog(mu_prop, sigma_prop, delta)

  if (feasible_prop) {
    BQ_prop <- gmm_BQr(x, pi_new, mu_prop, sigma_prop,
                       latent = latent, bw = bw, delta = delta)

    # feasible + non-decrease -> accept
    if (is.finite(BQ_old) && is.finite(BQ_prop) && (BQ_prop >= BQ_old)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
    # if old is not finite but new is finite, accept
    if (!is.finite(BQ_old) && is.finite(BQ_prop)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
  }

  # 2) backtracking: mu_try = (1-a)mu_old + a*mu_prop
  alpha <- alpha0
  mu_try <- mu_old
  sigma_try <- sigma_old
  BQ_try <- BQ_old

  for (bt in 1:max_backtrack) {
    alpha <- alpha * alpha_shrink

    mu_try <- lapply(seq_len(K), function(k) {
      (1 - alpha) * as.numeric(mu_old[[k]]) + alpha * as.numeric(mu_prop[[k]])
    })

    sigma_try <- gmm_sigma_update(x, mu_try, latent)

    # feasibility first
    if (!is_feasible_sumlog(mu_try, sigma_try, delta)) next

    BQ_try <- gmm_BQr(x, pi_new, mu_try, sigma_try,
                      latent = latent, bw = bw, delta = delta)

    if (is.finite(BQ_try) && is.finite(BQ_old) && (BQ_try >= BQ_old)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
    if (!is.finite(BQ_old) && is.finite(BQ_try)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
  }

  # 3) fallback: null move
  return(list(mu = mu_old, sigma = sigma_old, BQ = BQ_old, alpha = 0.0))
}


## ------------------------------------------------------------
## Helper func
## ------------------------------------------------------------

gmm_estep_annealed <- function(x, theta, r = 1) {
  x <- as.matrix(x)
  n <- nrow(x)
  K <- length(theta$pi)

  # log_resp[i,k] = r * (log pi_k + log N(x_i | mu_k, Sigma_k))
  log_resp <- matrix(NA_real_, nrow = n, ncol = K)

  for (k in 1:K) {
    log_comp <- dmvnorm_log(x, theta$mu[[k]], theta$Sigma[[k]])
    log_resp[, k] <- r * (log(theta$pi[k]) + log_comp)
  }

  # w_ik = exp(log_resp_ik)
  w <- exp(log_resp)

  # gamma_ik = w_ik / sum_j w_ij
  gamma <- w / rowSums(w)

  gamma
}

gmm_pi_update <- function(x, latent){
  N <- nrow(x)
  pi <- colSums(latent) / N
  as.numeric(pi)
}

gmm_mu_update <- function(x, latent){
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  K <- ncol(latent)

  Nk <- colSums(latent)          # K-vector
  mu <- vector("list", K)

  for(k in 1:K){
    w <- latent[, k]             # n-vector
    mu[[k]] <- as.numeric(colSums(x * w) / Nk[k])
  }

  mu
}

gmm_sigma_update <- function(x, mu, latent){
  x <- as.matrix(x)
  d <- ncol(x)
  K <- ncol(latent)

  Nk <- colSums(latent)
  Sigma <- vector("list", K)

  for(k in 1:K){
    xc <- sweep(x, 2, mu[[k]], "-")     # n x d
    w  <- latent[, k]                   # n-vector
    Sigma[[k]] <- (t(xc) %*% (xc * w)) / Nk[k]  # d x d
  }

  Sigma
}

gmm_em_at_r <- function(x, theta_init, r = 1,
                        tol = 1e-3, max_iter = 200, verbose = FALSE) {
  theta <- theta_init
  trace <- data.frame(iter=integer(0), r=numeric(0), Qr=numeric(0))

  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma

  Qr_old <- gmm_Qr(x, pi, mu, sigma, latent=gmm_estep_annealed(x, theta, r))

  # write trace
  trace <- rbind(trace, data.frame(iter=0, r=r, Qr=Qr_old))

  for (it in 1:max_iter) {
    # E-step
    latent = gmm_estep_annealed(x, theta, r)

    # M-step
    pi_new    = gmm_pi_update(x,latent)
    mu_new    = gmm_mu_update(x,latent)
    sigma_new = gmm_sigma_update(x,mu_new,latent)

    # Comupte new Qr
    Qr_new <- gmm_Qr(x, pi_new, mu_new, sigma_new, latent)

    # Update parameter
    theta$pi <- pi_new
    theta$mu <- mu_new
    theta$Sigma <- sigma_new

    # write trace
    trace <- rbind(trace, data.frame(iter=it, r=r, Qr=Qr_new))
    if (verbose) cat("iter:", it, "Qr:", Qr_new, "dQr:", Qr_new - Qr_old, "\n")

    # Break rule
    if(abs(Qr_old-Qr_new)<tol){
      break}else{Qr_old=Qr_new}
  }


  list(theta = theta, trace = trace)
}

gmm_em_at_r_bw <- function(x, theta,
                           r, bw,
                           delta,
                           tol = 1e-8, max_iter = 1000,
                           verbose = FALSE) {
  # optimization settings for inner mu solver
  maxiter_optim <- 200
  tol_optim <- 1e-4

  # trace of objective values
  trace <- data.frame(iter=integer(0), r=numeric(0), bw=numeric(0), Qr=numeric(0))

  # unpack parameters
  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma

  K <- length(mu)
  mu_new <- vector("list", K)

  ## initial barrier-augmented Q
  BQr_old <- gmm_BQr(
    x, pi, mu, sigma,
    latent = gmm_estep_annealed(x, theta, r),
    bw = bw, delta = delta
  )

  trace <- rbind(trace, data.frame(iter=0, r=r, bw=bw, Qr=BQr_old))

  for (it in 1:max_iter) {
    ## E-step
    latent <- gmm_estep_annealed(x, theta, r)

    ## M-step: update pi
    pi_new <- gmm_pi_update(x, latent)

    ## M-step: update mu by Gauss–Seidel with sumlog barrier solver
    mu_prop <- mu
    for (k in seq_len(K)) {
      mu_prop[[k]] <- solve_BQ_sumlog_mubarrier_newton(
        k = k, x = x,
        mu = mu_prop,
        sigma = sigma,
        latent = latent,
        delta = delta, bw = bw,
        maxiter = maxiter_optim, tol = tol_optim,
        alpha_shrink = 0.5
      )
    }

    ## accept update with feasibility and BQ non-decrease check
    acc <- accept_mu_sigma_by_damping_feasible_sumlog(
      mu_old = mu, sigma_old = sigma,
      x = x, pi_new = pi_new, latent = latent,
      mu_prop = mu_prop,
      bw = bw, delta = delta,
      BQ_old = BQr_old,
      alpha0 = 1.0, alpha_shrink = 0.5, max_backtrack = 20
    )

    mu_new    <- acc$mu
    sigma_new <- acc$sigma
    BQr_new   <- acc$BQ

    ## update parameters
    theta$pi    <- pi_new
    theta$mu    <- mu_new
    theta$Sigma <- sigma_new

    ## refresh locals
    pi    <- pi_new
    mu    <- mu_new
    sigma <- sigma_new

    ## update trace
    trace <- rbind(trace, data.frame(iter=it, r=r, bw=bw, Qr=BQr_new))

    if (verbose) {
      cat("iter:", it,
          "BQr(sumlog):", BQr_new,
          "dBQr:", BQr_new - BQr_old,
          "\n")
    }

    ## stop if improvement is small
    if (abs(BQr_old - BQr_new) < tol) {
      break
    } else {
      BQr_old <- BQr_new
    }
  }

  list(theta = theta, trace = trace)
}

make_trace_entry <- function(step, out, trace_full = TRUE, ...) {
  meta <- list(...)
  if (trace_full) {
    return(c(list(step = step), meta, list(out = out)))
  } else {
    return(c(list(
      step = step,
      BQ_final = tail(out$trace$Qr, 1),
      iters = nrow(out$trace) - 1
    ), meta))
  }
}

gmm_delta_DKL <- function(x, theta0, theta1, r) {
  x <- as.matrix(x)
  n <- nrow(x)

  # standard posterior (r=1) for theta0, theta1
  gamma0 <- gmm_estep_annealed(x, theta0, r = 1)
  gamma1 <- gmm_estep_annealed(x, theta1, r = 1)

  # annealed posterior weight (r=r) under theta0
  w <- gmm_estep_annealed(x, theta0, r = r)

  # per-i contribution: sum_k w_ik * (log gamma0_ik - log gamma1_ik)
  contrib <- rowSums(w * (log(gamma0) - log(gamma1)))

  # empirical mean over i
  mean(contrib)
}


gmm_find_adap_dhem_at_r <- function(x, theta, r = 1,
                                    bw = 1e-2,
                                    eta = 0.1,
                                    delta=0.1,
                                    tol = 1e-3, max_iter = 200,
                                    verbose = FALSE) {

  # optimization settings for inner mu solver
  maxiter_optim <- 200
  tol_optim <- 1e-4

  # trace of objective values
  trace <- data.frame(iter=integer(0), r=numeric(0), bw=numeric(0), Qr=numeric(0))

  # store initial theta for fallback
  theta0 = theta

  # unpack parameters
  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma

  K <- length(mu)
  mu_new <- vector("list", K)

  ## initial BQ at theta0
  BQr_old <- gmm_BQr(
    x, pi, mu, sigma,
    latent = gmm_estep_annealed(x, theta, r),
    bw = bw, delta = delta
  )

  trace <- rbind(trace, data.frame(iter=0, r=r, bw=bw, Qr=BQr_old))

  for (it in 1:max_iter) {

    ## ----- E-step -----
    latent <- gmm_estep_annealed(x, theta, r)

    ## ----- M-step -----
    pi_new <- gmm_pi_update(x, latent)

    ## ----- M-step: mu update (Gauss–Seidel with sumlog barrier solver) -----
    mu_prop <- mu
    for (k in seq_len(K)) {
      mu_prop[[k]] <- solve_BQ_sumlog_mubarrier_newton(
        k = k, x = x,
        mu = mu_prop,
        sigma = sigma,
        latent = latent,
        delta = delta, bw = bw,
        maxiter = maxiter_optim, tol = tol_optim,
        alpha_shrink = 0.5
      )
    }

    ## ----- Backtracking acceptance (feasible + BQ non-decrease) -----
    acc <- accept_mu_sigma_by_damping_feasible_sumlog(
      mu_old = mu, sigma_old = sigma,
      x = x, pi_new = pi_new, latent = latent,
      mu_prop = mu_prop,
      bw = bw, delta = delta,
      BQ_old = BQr_old,
      alpha0 = 1.0, alpha_shrink = 0.5, max_backtrack = 20
    )

    mu_new    <- acc$mu
    sigma_new <- acc$sigma
    BQr_new   <- acc$BQ

    ## update parameters
    theta$pi    <- pi_new
    theta$mu    <- mu_new
    theta$Sigma <- sigma_new

    ## refresh locals
    pi    <- pi_new
    mu    <- mu_new
    sigma <- sigma_new

    ## break if infeasible
    feasible_now <- is.finite(BQr_new)
    if(!is.finite(BQr_new)){break}

    ## update trace
    trace <- rbind(trace, data.frame(iter=it, r=r, bw=bw, Qr=BQr_new))

    ## ----- Acceptance condition (adaptive criterion) -----
    dkl_r <- gmm_delta_DKL(x, theta0, theta, r = r)
    dkl_1 <- gmm_delta_DKL(x, theta0, theta, r = 1)
    target <- eta * dkl_1
    acc_now <- (dkl_r >= target)

    trace <- rbind(trace,data.frame(iter = it, r = r, bw = bw,Qr = BQr_new))

    if (verbose) {
      cat("iter:", it,
          "BQ:", BQr_new,
          "dBQ:", BQr_new - BQr_old,
          "feasible:", feasible_now,
          "dkl_r:", dkl_r,
          "target:", target,
          "accept:", acc_now, "\n")
    }

    ## return if acceptance condition satisfied
    if (acc_now) {
      return(list(theta = theta,
                  trace = trace,
                  accept = TRUE,
                  bw = bw,
                  BQ_old = BQr_old,
                  BQ_new = BQr_new,
                  dkl_r = dkl_r,
                  dkl_1 = dkl_1,
                  target = target))
    }

    ## stop if improvement is small
    if (abs(BQr_new - BQr_old) < tol) break

    ## update for next iteration
    BQr_old <- BQr_new
  }

  ## return initial theta if no acceptable update found
  list(theta = theta0,
       trace = trace,
       accept = FALSE,
       bw = bw,
       BQ_old = BQr_old,
       BQ_new = BQr_new,
       dkl_r = dkl_r,
       dkl_1 = dkl_1,
       target = target)
}




## ------------------------------------------------------------
## Algorithms
## ------------------------------------------------------------
gmm_em <- function(x, theta_init,
                   tol = 1e-4, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  out = gmm_em_at_r(x, theta_init, r = 1, tol = tol, max_iter = max_iter, verbose = verbose)

  trace <- make_trace_entry(step = 1, out = out, trace_full = trace_full, r = 1)

  list(theta = out$theta, trace = trace)
}

gmm_daem <- function(x, theta_init,
                     r_init = 0.2, n_steps = 50,
                     tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  theta <- theta_init
  trace <- vector("list", n_steps)
  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))

  for (s in seq_along(r_grid)) {
    r <- r_grid[s]
    out <- gmm_em_at_r(x, theta, r = r, tol = tol, max_iter = max_iter, verbose = verbose)
    theta <- out$theta

    # r by write
    trace[[s]] <- make_trace_entry(step = s, out = out, trace_full = trace_full, r = r)
  }
  list( theta = theta,trace = trace)
}


gmm_dhem <- function(x, theta_init,
                     r_init = 0.2, n_steps = 50,
                     bw_init = 1, bw_end = 1e-10,
                     delta = 0.1,
                     tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {

  theta <- theta_init
  trace <- vector("list", n_steps)

  r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))                 # r -> 1
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))      # bw -> bw_end

  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]
    bw <- bw_grid[s]

    out <- gmm_em_at_r_bw(
      x, theta,
      r = r, bw = bw, delta = delta,
      tol = tol, max_iter = max_iter, verbose = verbose
    )

    theta <- out$theta

    trace[[s]] <- make_trace_entry(step = s, out = out, trace_full = trace_full, r = r, bw = bw)
  }

  list(theta = theta, trace = trace  )
}

gmm_barrier <- function(x, theta_init,
                        n_steps = 50,
                        # r_init = 0.2,
                        bw_init = 1, bw_end = 1e-10,
                        delta = 0.1,
                        tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {

  theta <- theta_init
  trace <- vector("list", n_steps)

  # r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))

  r = 1
  for (s in seq_len(n_steps)) {
    # r  <- r_grid[s]
    bw <- bw_grid[s]

    out <- gmm_em_at_r_bw(
      x, theta,
      r = r, bw = bw, delta = delta,
      tol = tol, max_iter = max_iter, verbose = verbose
    )

    theta <- out$theta

    trace[[s]] <- make_trace_entry(step=s, out=out, trace_full=trace_full, r=1, bw=bw)
  }

  list(theta = theta, trace = trace  )
}



gmm_adaptive_dhem <- function(x, theta,
                              r_init = 0.2, n_steps = 50,
                              bw_init = 1e-2,
                              delta = 0.1,
                              eta = 0.1,
                              tol_BQ = 1e-4,
                              tol_find = 1e-3,
                              max_iter_inner = 50,
                              max_iter_find = 200,
                              verbose = FALSE,
                              trace_full = TRUE) {

  x <- as.matrix(x)

  # geometric schedule for r
  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))

  # barrier weight (monotonically decreasing)
  bw <- bw_init

  trace <- vector("list", n_steps)

  for (s in seq_along(r_grid)) {

    r <- r_grid[s]

    if (verbose) cat("\n===== step", s, " r =", r, " bw =", bw, "=====\n")

    BQ_prev_acc <- -Inf
    out_last <- NULL
    acc_count <- 0L

    # inner loop at fixed r
    for (inner in seq_len(max_iter_inner)) {

      out <- gmm_find_adap_dhem_at_r(
        x = x,
        theta = theta,
        r = r,
        bw = bw,
        eta = eta,
        delta = delta,
        tol = tol_find,
        max_iter = max_iter_find,
        verbose = verbose
      )
      out_last <- out

      # accept1 failure → move to next r
      if (!out$accept) {
        if (verbose) cat("  inner:", inner, " accept1=FALSE -> next r\n")
        break
      }

      theta_new <- out$theta
      BQ_new    <- out$BQ_new

      # check barrier condition (2)
      cap <- cap_bw_by_barrier(
        theta0 = theta,
        theta1 = theta_new,
        eta = eta,
        dkl_1 = out$dkl_1,
        bw = bw,
        delta = delta
      )

      # shrink bw if violated and retry at same r
      if (isTRUE(cap$violated)) {
        if (verbose) {
          cat("  inner:", inner,
              " barrier violated -> bw:",
              bw, "->", cap$bw_new, "\n")
        }
        bw <- cap$bw_new
        next
      }

      # accept update
      theta <- theta_new
      acc_count <- acc_count + 1L

      # stop if BQ improvement is small
      if (is.finite(BQ_new) && abs(BQ_new - BQ_prev_acc) < tol_BQ) {
        if (verbose) cat("  inner:", inner, " small ΔBQ -> break\n")
        BQ_prev_acc <- BQ_new
        break
      }
      BQ_prev_acc <- BQ_new
    }

    # store trace
    if (trace_full) {
      trace[[s]] <- list(step = s, r = r, bw = bw,
                         accepted_moves = acc_count,
                         out_last = out_last,
                         theta = theta)
    } else {
      if (is.null(out_last)) {
        trace[[s]] <- list(step = s, r = r, bw = bw,
                           accepted_moves = 0L,
                           accept_last = FALSE,
                           BQ_final = NA_real_,
                           iters_last = 0L)
      } else {
        BQ_final <- tail(out_last$trace$BQ, 1)
        trace[[s]] <- list(step = s, r = r, bw = bw,
                           accepted_moves = acc_count,
                           accept_last = out_last$accept,
                           BQ_final = BQ_final,
                           iters_last = nrow(out_last$trace) - 1L)
      }
    }
  }

  list(theta = theta, bw_final = bw, trace = trace)
}



## -----------------------------
## distances + matching (minimal)
## -----------------------------
l2_dist <- function(a, b) sqrt(sum((as.numeric(a) - as.numeric(b))^2))
fro_dist <- function(A, B) sqrt(sum((A - B)^2))

match_by_mu <- function(theta_hat, theta_true) {
  K <- length(theta_true$pi)
  C <- matrix(0, K, K)
  for (k in 1:K) for (j in 1:K) {
    C[k, j] <- l2_dist(theta_hat$mu[[k]], theta_true$mu[[j]])
  }
  used_hat <- rep(FALSE, K)
  perm_hat_to_true <- integer(K)
  for (j in 1:K) {
    cand <- which(!used_hat)
    kbest <- cand[which.min(C[cand, j])]
    perm_hat_to_true[j] <- kbest
    used_hat[kbest] <- TRUE
  }
  th <- list(
    pi = theta_hat$pi[perm_hat_to_true],
    mu = theta_hat$mu[perm_hat_to_true],
    Sigma = theta_hat$Sigma[perm_hat_to_true]
  )
  list(theta = th, perm = perm_hat_to_true)
}

compute_err3 <- function(theta_hat, theta_true) {
  K <- length(theta_true$pi)
  th <- match_by_mu(theta_hat, theta_true)$theta

  pi_abs <- mean(abs(th$pi - theta_true$pi))

  mu_L2 <- mean(vapply(seq_len(K),
                       function(j) l2_dist(th$mu[[j]], theta_true$mu[[j]]),
                       numeric(1)))

  Sigma_F <- mean(vapply(seq_len(K),
                         function(j) fro_dist(th$Sigma[[j]], theta_true$Sigma[[j]]),
                         numeric(1)))

  c(pi = pi_abs, mu = mu_L2, Sigma = Sigma_F)
}


gmm_bw_init_mu <- function(x, theta_init, r_init = 1, delta = 0.1,
                           tau = 0.05, eps = 1e-12,
                           norm_type = c("l2", "linf")) {
  norm_type <- match.arg(norm_type)

  x <- as.matrix(x)
  theta <- theta_init
  latent <- gmm_estep_annealed(x, theta, r = r_init)

  K <- length(theta$mu)
  vals <- numeric(K)

  vec_norm <- function(v) {
    if (norm_type == "l2") sqrt(sum(v^2)) else max(abs(v))
  }

  for (k in seq_len(K)) {
    w  <- latent[, k]
    Nk <- sum(w)
    sk <- colSums(x * w)
    mk <- as.numeric(theta$mu[[k]])

    gQk <- sk - Nk * mk

    mu_others <- theta$mu[-k]
    uk <- Reduce(`+`, lapply(mu_others, function(ml) mk - as.numeric(ml)))

    Sig_k <- theta$Sigma[[k]]
    Rk <- chol(Sig_k)

    Dk <- 0
    for (ml in mu_others) {
      diff <- mk - as.numeric(ml)
      y <- backsolve(Rk, diff, transpose = TRUE)
      Dk <- Dk + sum(y^2)
    }
    Dk <- Dk - delta

    if (!is.finite(Dk) || Dk <= 0) {
      vals[k] <- 0
      next
    }

    vals[k] <- (tau / 2) * vec_norm(gQk) * Dk / (vec_norm(uk) + eps)
  }

  min(vals)
}

## -----------------------------
## one simulation summary (requested format)
## -----------------------------

one_sim_summary <- function(simulation_num = 1,
                            n, K, theta_true,
                            tol_limit = 1e-6,
                            r_init = 0.2, n_steps = 20,
                            bw_init = 1e-2, bw_end = 1e-10,
                            delta = 0.5,
                            max_iter = 1000,eta=0.1,
                            verbose = FALSE) {

  # random seed for this simulation
  seed <- sample.int(1e9, 1)

  # generate data and initial parameters
  tmp <- make_data_and_init(n, K, theta_true, delta)
  x <- tmp$x
  theta0 <- tmp$theta_init


  # barrier initial (mu)
  bw_init_mu <- gmm_bw_init_mu(
    x = x,
    theta_init = theta0,
    r_init = r_init,
    delta = delta,
    tau = 0.01,        # 핵심 tuning
    norm_type = "l2"
  )

  bw_init = bw_init_mu

  # safe wrapper: return NULL if error occurs
  safe_run <- function(expr) tryCatch(expr, error = function(e) NULL)

  # run all methods
  res <- list(
    EM = safe_run(gmm_em(x, theta0, tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    DAEM = safe_run(gmm_daem(x, theta0,
                             r_init = r_init, n_steps = n_steps,
                             tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    DHEM = safe_run(gmm_dhem(x, theta0,
                             r_init = r_init, n_steps = n_steps,
                             bw_init = bw_init, bw_end = bw_end,
                             delta = delta,
                             tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    Barrier = safe_run(gmm_barrier(x, theta0,
                                   n_steps = n_steps,
                                   bw_init = bw_init, bw_end = bw_end,
                                   delta = delta,
                                   tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    AdaptiveDHEM = safe_run(gmm_adaptive_dhem(x, theta0,
                                              r_init = r_init, n_steps = n_steps,
                                              bw_init = bw_init,
                                              delta = delta,
                                              eta = 0.1,
                                              tol_BQ = tol_limit, tol_find = 1e-3,
                                              max_iter_inner = max_iter, max_iter_find = 200,
                                              verbose = verbose, trace_full = TRUE))
  )

  # extract estimated parameters
  theta_list <- lapply(res, function(o) if (is.null(o)) NULL else o$theta)

  methods <- names(theta_list)
  paras <- c("pi", "mu", "Sigma")

  rows <- list()
  idx <- 0L

  for (m in methods) {
    th_hat <- theta_list[[m]]
    succ <- if (is.null(th_hat)) 0L else 1L

    # if failed, record NA errors
    if (succ == 0L) {
      for (p in paras) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          simulation_num = simulation_num,
          seed = seed,
          method = m,
          para = p,
          err = NA_real_,
          succ = succ,
          stringsAsFactors = FALSE
        )
      }
    } else {
      # compute estimation error for each parameter
      e <- compute_err3(th_hat, theta_true)
      for (p in paras) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          simulation_num = simulation_num,
          seed = seed,
          method = m,
          para = p,
          err = as.numeric(e[p]),
          succ = succ,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # combine results into a single data frame
  do.call(rbind, rows)
}


## ----------------------------------------------------------------
## Section 2: Weibull Mixture Model (WMM) functions
## Source: wmm_functions.R
## ----------------------------------------------------------------
library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(patchwork)
library(purrr)

############################
# Help functions
############################

diffB_onlyB = function(beta,event_vec,time_vec,latentZ_mat,j){
  sum(latentZ_mat[,j]*event_vec)/beta +
    sum(latentZ_mat[,j]*event_vec*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)*sum(latentZ_mat[,j]*(time_vec^beta)*log(time_vec))/sum(latentZ_mat[,j]*(time_vec^beta))
}

barrierFunc_1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1)+(1/beta -1/(1-beta))*(bw)
  return(result)
}

barrierFunc_3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3)+bw*(1/(beta-1))
  return(result)
}

barrier_beta1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1))

  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps = 1e-12
  if(bw==0){
    while(!isna_diffbeta(maxRange)){
      maxRange=maxRange + min(maxRange*1.01,10)
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=1)*
      diffB_onlyB(1e-12,event_vec,time_vec,latentZ_mat, j=1)<0) break
    }
      result <- uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1),
      interval = c(0,maxRange),tol=1e-10)
      return(result$root)
  }

  result <- uniroot(function(beta) barrierFunc_1(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(eps,1-eps),tol=1e-10)
  return(result$root)
}

barrier_beta3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3))

  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps= 1e-3
  if(bw==0){
    while(!isna_diffbeta(maxRange)&&!isna_diffbeta(eps)){
      maxRange=maxRange*1.01
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
      diffB_onlyB(eps,event_vec,time_vec,latentZ_mat, j=3)<0) break
    }
    result = uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3),
    interval = c(eps, maxRange),tol=1e-10)
    return(result$root)
  }

  while(!isna_diffbeta(maxRange)){
    maxRange=maxRange + min(maxRange*1.01,10)
    if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
    diffB_onlyB(1,event_vec,time_vec,latentZ_mat, j=3)<0) break
  }
  result = uniroot(function(beta) barrierFunc_3(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}

barrier_safe_wrapper1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta1(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

barrier_safe_wrapper3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta3(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

weibull_estep_annealed <- function(df, pi, lambda, beta, r = 1) {
  t     <- as.numeric(df$time)
  event <- as.numeric(df$event)

  n <- length(t)
  K <- length(pi)

  logt <- log(t)
  log_resp <- matrix(NA_real_, nrow = n, ncol = K)

  for (k in 1:K) {
    # Weibull component k:
    # f_k(t) = lambda_k * beta_k * t^(beta_k-1) * exp(-lambda_k * t^beta_k)
    # S_k(t) = exp(-lambda_k * t^beta_k)
    #
    # log L_ik = event_i * log f_k(t_i) + (1-event_i) * log S_k(t_i)
    #         = event_i*(log lambda_k + log beta_k + (beta_k-1)log t_i) - lambda_k t_i^beta_k
    logLik_ik <- event * (log(lambda[k]) + log(beta[k]) + (beta[k] - 1) * logt) -
      lambda[k] * (t ^ beta[k])

    # annealed responsibilities: gamma_ik ∝ exp( r * (log pi_k + logLik_ik) )
    log_resp[, k] <- r * (log(pi[k]) + logLik_ik)
  }

  # numeric stabilization and normalization across k
  rowmax <- apply(log_resp, 1, max)
  w <- exp(log_resp - rowmax)
  gamma <- w / rowSums(w)

  gamma
}

wmm_delta_DKL <- function(df, theta0, theta1, r) {

  # standard posterior (r=1) for theta0, theta1
  gamma0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  gamma1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)

  # annealed posterior weight (r=r) under theta0
  w <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = r)

  # per-i contribution: sum_k w_ik * (log gamma0_ik - log gamma1_ik)
  contrib <- rowSums(w * (log(gamma0) - log(gamma1)))

  mean(contrib)
}

wmm_DKL <- function(df, theta0, theta1) {
  g0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  g1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)
  mean(rowSums(g0 * (log(g0) - log(g1))))
}

wmm_bar_value <- function(beta) {
  log(beta[1]) + log(1 - beta[1]) + log(beta[3] - 1)
}

wmm_bar_diff <- function(theta_new, theta_old) {
  wmm_bar_value(theta_new$beta) - wmm_bar_value(theta_old$beta)
}

wmm_lambda_init <-function(time_vec,event_vec,beta_vec,ratio1,ratio3){
  surv_obj <- Surv(time = time_vec, event = event_vec)
  fit <- survfit(surv_obj ~ 1)

  times <- fit$time
  surv_probs <- fit$surv

  cumhaz <- -log(surv_probs)

  delta_time <- diff(c(0, times))
  delta_hazard <- diff(c(0, cumhaz))
  hazard_rate <- delta_hazard / delta_time

  df_haz = data.frame(time=unique(time_vec),hazard=hazard_rate,
                      time1= beta_vec[1] * unique(time_vec)^(beta_vec[1] - 1),
                      time3= beta_vec[3] * unique(time_vec)^(beta_vec[3] - 1)
                      )
  df_haz[(nrow(df_haz)),"hazard"] = df_haz[(nrow(df_haz)-1),"hazard"]
  n = nrow(df_haz)
  subset_df1 = df_haz %>% filter(time<max(times)*ratio1)
  subset_df3 = df_haz %>% filter(time>max(times)*ratio3)
  subset_df2 = df_haz %>% filter(time>max(times)*ratio1 & time<max(times)*ratio3 )
  fit1 <- lm(hazard ~ 0 + time1, data = subset_df1)  # '0 +'는 intercept 제외
  fit3 <- lm(hazard ~ 0 + time3, data = subset_df3)  # '0 +'는 intercept 제외

  # 결과 확인
  lambda_est1 <- coef(fit1)[1]
  lambda_est3 <- coef(fit3)[1]
  lambda_est2 = mean(subset_df2[,"hazard"])
  lambda_vec = c(lambda_est1,lambda_est2,lambda_est3)
  return(lambda_vec)
}


############################
# Algorithm
############################

wmm_EM <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   tol = 1e-6,verbose=FALSE
                  ){

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", length(maxGEMiter))

  bw = 0
  for (it in 1:maxGEMiter) {
    ### M-step ###
    new_pi = colSums(gamma)/N

    new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
    new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = 1
    new_beta = c(new_beta1,new_beta2,new_beta3)

    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))

    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    beta=new_beta; pi=new_pi; lambda=new_lambda;

    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
    ### Stopping rule ###
    if(parameter_diff<tol || it==maxGEMiter){
      if(verbose) cat("EM ","[",it,"]"," beta :",beta ,"\n")
      break
    }
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=1)

    trace[[it]] <- list(pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
  }

  return(list(trace = trace, lambda = lambda, beta = beta))
}


wmm_DAEM <- function(df,
                     theta,
                     maxGEMiter = 1e+3,
                     nsteps = 50,
                     r_init = 0.1,
                     r_end = 1,
                     tol=1e-6,verbose=FALSE
                    ){
  method = "DAEM"
  r_grid <- exp(seq(log(r_init), log(r_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  bw = 0
  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]

    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N

      new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)

      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))

      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){

        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }

  return(list(trace = trace, pi = pi, lambda = lambda, beta = beta))
}

wmm_BM <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   nsteps = 50,
                   bw_init = 1e-1,
                   bw_end = 1e-5,
                   tol=1e-6,verbose=FALSE
                  ){
  method="BM"
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  r=1

  for( hyperIter in 1:nsteps){

    bw=bw_grid[hyperIter]

    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N

      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)

      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))

      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){

        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM = function(df,
                    theta,
                    maxGEMiter=1e+3,
                    nsteps=1e+2,
                    r_init=0.1,
                    r_end=1,
                    bw_init=1e-1,
                    bw_end=1e-5,
                    tol=1e-6,verbose=FALSE
  ){
  method = "DHEM"
  r_grid  <- exp(seq(log(r_init), log(r_end), length.out = nsteps))
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]
    bw=bw_grid[hyperIter]

    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N

      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)

      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))

      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2+(pi-new_pi)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){

        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM_adaptive <- function(df,
                              theta,
                              maxGEMiter=1e+3,
                              nsteps=1e+2,
                              r_init=0.1,
                              r_end=1,
                              bw_init=1e-1,
                              eta=0.1,
                              tol=1e-6,verbose=FALSE
                            ){
  method = "adapDHEM"
  # r schedule (r -> 1)
  r_grid <- exp(seq(log(r_init), 0, length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  bw = bw_init
  for( hyperIter in 1:nsteps){
    r <- r_grid[hyperIter]
    #bw = bw*0.1
    for( gemIter in 1:maxGEMiter){
    ### M-step ###
    new_pi = colSums(gamma)/N

    new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
    new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = 1
    new_beta = c(new_beta1,new_beta2,new_beta3)

    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))

    theta0 = list(pi=pi,beta=beta,lambda=lambda)
    theta1 = list(pi=new_pi,beta=new_beta,lambda=new_lambda)

    # Checking ACC
    acc1 = FALSE;acc2 = FALSE;
    deltaDKL = wmm_delta_DKL(df,theta0,theta1,r)
    DKL      = wmm_DKL(df,theta0,theta1)
    deltaB   = wmm_bar_diff(theta1,theta0)
    #cat(deltaDKL-bw*wmm_bar_diff(theta1,theta0),"\n")
    if(deltaDKL-bw*deltaB<0){
      # Acc 1st test
      if(!is.finite(deltaDKL)||!is.finite(DKL)||DKL<0) break
      if(deltaDKL<eta*DKL) {
        #cat(deltaDKL,eta*DKL,"\n")
        break}else{acc1=TRUE}
      # Acc 2nd test
      # if(bw*abs(deltaB)>eta*DKL){
      if(bw*abs(deltaB)>0.5*eta*DKL){
        bw = min(bw,0.5*eta*DKL/abs(deltaB))
        next
      }else{acc2 = TRUE}
    }else{
      acc1=TRUE;acc2=TRUE;
    }
    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    #print(parameter_diff)
    beta=new_beta; pi=new_pi; lambda=new_lambda;
    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)

    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    ### Stopping rule ###
    if(parameter_diff<tol || gemIter==maxGEMiter){
      if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3," para diff: ", parameter_diff,"\n")
      break
    }
    }

    ## Record this hyperIter's outcome only once, using the status of the
    ## LAST gemIter attempt actually made at this r. If that final attempt
    ## was rejected (acc1 or acc2 FALSE) -- even after earlier gemIter steps
    ## within this same hyperIter were accepted -- the whole hyperIter is
    ## discarded and last_acc/trace keep the previous hyperIter's state.
    ## This matches the original (pre-refactor) implementation and is what
    ## reproduces the paper's Adaptive DHEM numbers: a hyperIter whose final
    ## attempt fails typically means the barrier has been driven towards the
    ## boundary faster than genuine likelihood improvement justifies, so its
    ## (possibly boundary-hugging) intermediate progress should not be kept
    ## as the reported estimate.
    if(acc1&&acc2){
      trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
      last_acc = list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
    }
  }

  return(list(trace = trace,pi = last_acc$pi, lambda = last_acc$lambda, beta = last_acc$beta,bw=last_acc$bw,dQbeta1 = last_acc$dQbeta1, dQbeta3=last_acc$dQbeta3))
}

make_result_df <- function(fit, method, K = 3) {
  non_null_idx  <- which(!sapply(fit$trace, is.null))
  valid_traces  <- fit$trace[non_null_idx]

  data.frame(
    method  = method,
    r       = sapply(valid_traces, function(x) if (!is.null(x$r))  x$r  else 1),
    bw      = sapply(valid_traces, function(x) if (!is.null(x$bw)) x$bw else 0),
    pi1     = sapply(valid_traces, function(x) x$pi[1]),
    pi2     = sapply(valid_traces, function(x) x$pi[2]),
    pi3     = sapply(valid_traces, function(x) x$pi[3]),
    beta1   = sapply(valid_traces, function(x) x$beta[1]),
    beta2   = sapply(valid_traces, function(x) x$beta[2]),
    beta3   = sapply(valid_traces, function(x) x$beta[3]),
    lambda1 = sapply(valid_traces, function(x) x$lambda[1]),
    lambda2 = sapply(valid_traces, function(x) x$lambda[2]),
    lambda3 = sapply(valid_traces, function(x) x$lambda[3]),
    dQbeta1 = sapply(valid_traces, function(x) if (length(x$dQbeta1) == 1) x$dQbeta1 else NA_real_),
    dQbeta3 = sapply(valid_traces, function(x) if (length(x$dQbeta3) == 1) x$dQbeta3 else NA_real_)
  )
}

plot_beta_trace <- function(df, title = NULL,
                            axis_text_y_size  = 14,
                            axis_title_y_size = 14,
                            axis_text_x_size  = 14,
                            axis_title_x_size = 14,
                            title_size        = 20
                          ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,
      angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = r, y = beta1)) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1), xlim = c(r_init, 1)) +
    labs(x = "Annealing parameter", y = expression(beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = r, y = beta3)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1)) +
    labs(x = "Annealing parameter", y = expression(beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


plot_dQbeta_trace <- function(df, title = NULL,
                              axis_text_y_size  = 14,
                              axis_title_y_size = 14,
                              axis_text_x_size  = 14,
                              axis_title_x_size = 14,
                              title_size        = 20
                            ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle =0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = r, y = dQbeta1)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1))+
    labs(x = "Annealing parameter", y = expression(nabla~Q~beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = r, y = dQbeta3)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1))+
    labs(x = "Annealing parameter", y = expression(nabla~Q~beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


best_row <- function(df) {
  idx <- which.min(abs(df$dQbeta1) + abs(df$dQbeta3))
  df[idx, c("method","pi1","pi2","pi3","beta1","beta3","lambda1","lambda2","lambda3","dQbeta1","dQbeta3")]
}


plot_beta_trace_barrier <- function(df, title = NULL,
                            axis_text_y_size  = 14,
                            axis_title_y_size = 14,
                            axis_text_x_size  = 10,
                            axis_title_x_size = 14,
                            title_size        = 20
                          ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,
      angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  bw_upper <- max(1, max(df$bw, na.rm = TRUE))

  p1 <- ggplot(df, aes(x = bw, y = beta1)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, bw_upper)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, bw_upper)
)+
    labs(x = "Barrier parameter", y = expression(beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = bw, y = beta3)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, bw_upper)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, bw_upper)
)+
    labs(x = "Barrier parameter", y = expression(beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


plot_dQbeta_trace_barrier <- function(df, title = NULL,
                              axis_text_y_size  = 14,
                              axis_title_y_size = 14,
                              axis_text_x_size  = 10,
                              axis_title_x_size = 14,
                              title_size        = 20
                            ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle =0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  bw_upper <- max(1, max(df$bw, na.rm = TRUE))

  p1 <- ggplot(df, aes(x = bw, y = dQbeta1)) +
    geom_line() +
    # scale_x_reverse() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, bw_upper)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, bw_upper)
)+
    labs(x = "Barrier parameter", y = expression(nabla~Q~beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = bw, y = dQbeta3)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, bw_upper)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, bw_upper)
)+
    labs(x = "Barrier parameter", y = expression(nabla~Q~beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}

wmm_bw_init_beta <- function(df, pi_init, lambda_init, beta_init, r_init,
                             tau = 0.1) {
  z_annealed <- weibull_estep_annealed(
    df, pi_init, lambda_init, beta_init, r_init
  )

  g1 <- diffB_onlyB(beta_init[1], df$event, df$time, z_annealed, j = 1)
  g3 <- diffB_onlyB(beta_init[3], df$event, df$time, z_annealed, j = 3)

  bw1 <- tau * abs(g1) * min(beta_init[1], 1 - beta_init[1])
  bw3 <- tau * abs(g3) * (beta_init[3] - 1)

  min(bw1, bw3)
}


## ----------------------------------------------------------------
## Section 3: Zero-Inflated Poisson (ZIP) model functions
## Source: ZIP_functions.R
## ----------------------------------------------------------------
# ZIP data generator ------------------------------------------------------

rzip <- function(n, pi, lambda) {
  # Z=0 (structural zero) with prob pi
  # Z=1 (Poisson) with prob 1-pi
  z <- rbinom(n, size = 1, prob = 1 - pi)  # z=1 => Poisson component
  x <- z * rpois(n, lambda = lambda)       # if z=0 => 0
  x
}

compute_lambda_init <- function(y, set_pi_init) {
  w <- ifelse(y == 0, set_pi_init, 1)
  sum(w * y) / sum(w)
}

rel_err <- function(est, true) {
  abs(est - true) / true
}

abs_err <- function(est, true) {
  abs(est - true)
}

err <- function(est, true) {
  (est - true)
}

ratio_err <- function(est, true) {
  max(abs(est/true),1/abs(est/true))
}

log_err <- function(est, true) {
  abs(log(est)-log(true))
}


# Annealed ZIP E-step ------------------------------------------------------

zip_estep_annealed <- function(x, pi, lambda, r=1) {
  if (r <= 0) stop("r must be > 0")

  gamma <- numeric(length(x))

  # x > 0 : 반드시 Poisson 성분
  pos <- (x > 0L)
  gamma[pos] <- 1.0

  # x = 0 : annealed posterior
  z0 <- !pos
  if (any(z0)) {
    w1 <- (1 - pi) * exp(-lambda)  # Poisson에서 0이 나올 joint weight
    w0 <- pi                       # structural zero weight

    w1r <- w1^r
    w0r <- w0^r
    gamma[z0] <- w1r / (w0r + w1r)
  }
  # γ=(1−π)e−λ/(π+(1−π)e−λ)
  gamma
}

zip_mstep <- function(x, gamma) {
  n <- length(x)
  A <- sum(1 - gamma)  # structural-zero 기대 개수
  B <- sum(gamma)      # Poisson component 기대 개수

  pi_new <- A / n
  lambda_new <- sum(gamma * x) / B

  list(pi = pi_new, lambda = lambda_new)
}


zip_em <- function(x,
                   pi_init = 0.5,
                   lambda_init = NULL,
                   tol = 1e-10,
                   max_iter = 200
                  ) {
  n <- length(x)

  # 초기값
  pi <- pi_init
  if (is.null(lambda_init)) {
    lambda <- compute_lambda_init(x, set_pi_init = pi_init)
  } else {
    lambda <- lambda_init
  }

  trace <- data.frame(
    iter = integer(0),
    pi = numeric(0),
    lambda = numeric(0),
    d_pi = numeric(0),
    d_lambda = numeric(0)
  )

  for (it in seq_len(max_iter)) {
    # E-step (r=1)
    gamma <- zip_estep_annealed(x, pi = pi, lambda = lambda, r = 1)

    # M-step (사용자 제공 함수 그대로)
    th_new <- zip_mstep(x, gamma)
    pi_new <- th_new$pi
    lambda_new <- th_new$lambda

    # 변화량 기록
    d_pi <- abs(pi_new - pi)
    d_lam <- abs(lambda_new - lambda)

    trace[it, ] <- list(it, pi_new, lambda_new, d_pi, d_lam)

    # 수렴 체크
    if (max(d_pi, d_lam) < tol) {
      pi <- pi_new
      lambda <- lambda_new
      break
    }

    # 업데이트
    pi <- pi_new
    lambda <- lambda_new
  }

  list(
    theta = list(pi = pi, lambda = lambda),
    trace = trace
  )
}
# ============================================================
# 0) 기존 함수들이 이미 정의되어 있다고 가정:
#    - zip_estep_annealed(x, pi, lambda, r)
#    - zip_mstep(x, gamma)
#    - zip_dkl_std_post(x, theta0, theta1)
#    - zip_delta_dkl(x, theta0, theta1, r)
# ============================================================


# ============================================================
# 1) fixed-r에서 EM을 tol까지 돌리는 "wrapper" (새 유틸)
#    -> 기존 zip_estep_annealed + zip_mstep만 사용
# ============================================================
zip_em_at_r <- function(x, theta_init, r,
                        tol = 1e-10,
                        max_iter = 200) {
  theta <- theta_init

  trace <- data.frame(
    iter   = integer(0),
    r      = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )

  for (it in 1:max_iter) {
    gamma  <- zip_estep_annealed(x, pi = theta$pi, lambda = theta$lambda, r = r)
    theta1 <- zip_mstep(x, gamma)

    d_pi  <- abs(theta1$pi - theta$pi)
    d_lam <- abs(theta1$lambda - theta$lambda)

    theta <- theta1

    trace <- rbind(trace, data.frame(
      iter   = it,
      r      = r,
      pi     = theta$pi,
      lambda = theta$lambda,
      d_pi   = d_pi,
      d_lam  = d_lam
    ))

    if (max(d_pi, d_lam) < tol) break
  }

  list(theta = theta, trace = trace)
}


# ============================================================
# 2) DAEM: r-grid (log schedule) 각 r에서 zip_em_at_r로 수렴
#    -> 기존 zip_daem의 "1회 업데이트" 문제 해결
# ============================================================
zip_daem <- function(x,
                     theta_init,
                     r_init  = 0.2,
                     n_steps = 50,
                     tol = 1e-10,
                     max_iter = 200) {

  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))  # ends at 1

  theta <- theta_init

  trace <- data.frame(
    step   = integer(0),
    r      = numeric(0),
    iter   = integer(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )

  for (s in 1:n_steps) {
    r <- r_grid[s]

    res  <- zip_em_at_r(x, theta_init = theta, r = r, tol = tol, max_iter = max_iter)
    theta <- res$theta

    tr <- res$trace
    tr$step <- s
    trace <- rbind(trace, tr[, c("step","r","iter","pi","lambda","d_pi","d_lam")])
  }

  list(theta = theta, trace = trace, r_grid = r_grid)
}


# ============================================================
# 3) Adaptive DAEM:
#    - r <= r_switch: fixed phase (각 r에서 inner EM 수렴, 항상 accept)
#    - r >  r_switch: 후보 r들을 증가시키며 (rescan)
#         theta1(r) = inner EM 수렴 결과로 만들고
#         zip_delta_dkl >= eta * zip_dkl_std_post 이면 accept
#         아니면 다음 r로 넘어감
# ============================================================

zip_BQ <- function(x, theta, gamma, bw = 0, pi_min = 0) {
  pi <- theta$pi
  lam <- theta$lambda
  if (!is.finite(pi) || !is.finite(lam)) return(-Inf)
  if (pi <= 0 || pi >= 1 || lam <= 0) return(-Inf)
  if (bw > 0 && pi <= pi_min) return(-Inf)

  Q <- sum((1 - gamma) * log(pi) + gamma * log(1 - pi) +
             gamma * (x * log(lam) - lam))
  B <- if (bw > 0) bw * log(pi - pi_min) else 0
  Q + B
}


# --- Barrier pi-update (M-step subroutine) -------------------------------

zip_mstep_pi_barrier <- function(A, B, bw, pi_min) {
  # maximize: f(pi) = A log(pi) + B log(1-pi) + bw log(pi - pi_min)
  # domain: pi in (pi_min, 1)
  # FOC -> quadratic:
  # n*pi^2 - (n*(1+pi_min)+bw)*pi + A*pi_min = 0
  n <- A + B

  b <- n * (1 + pi_min) + bw
  c <- A * pi_min

  disc <- b*b - 4*n*c
  if (disc < 0) disc <- 0  # numeric guard

  r1 <- (b - sqrt(disc)) / (2*n)
  r2 <- (b + sqrt(disc)) / (2*n)

  eps <- 1e-12
  lo  <- pi_min + eps
  hi  <- 1 - eps

  # candidate set: roots that lie in (pi_min, 1)
  cand <- c(r1, r2)
  cand <- cand[is.finite(cand) & (cand > lo) & (cand < hi)]

  # if no root is admissible (rare), fall back to interior projection of A/n
  if (length(cand) == 0) {
    p0 <- A / n
    return(min(max(p0, lo), hi))
  }

  # evaluate objective and pick best root
  f <- function(pi) A*log(pi) + B*log(1 - pi) + bw*log(pi - pi_min)
  vals <- sapply(cand, f)
  cand[which.max(vals)]
}


# --- ZIP M-step with optional barrier on pi ------------------------------

zip_mstep_barrier <- function(x, gamma, use_barrier = FALSE, bw = 0, pi_min = 0.005) {
  n <- length(x)
  A <- sum(1 - gamma)
  B <- sum(gamma)

  # lambda update (same as usual)
  lambda_new <- sum(gamma * x) / sum(gamma)

  # pi update
  if (!use_barrier) {
    pi_new <- A / n
  } else {
    pi_new <- zip_mstep_pi_barrier(A, B, bw = bw, pi_min = pi_min)
  }

  list(pi = pi_new, lambda = lambda_new)
}

# ZIP barrier-EM -----------------------------------------------------------
# - r is fixed (default r=1 => standard posterior in E-step)
# - barrier schedule bw_init -> bw_end over n_steps
# - uses your existing: zip_estep_annealed(), zip_mstep_barrier()

zip_barrier <- function(x,
                        theta_init,
                        bw_init = 5,
                        bw_end  = 0.05,
                        pi_min  = 0.005,
                        n_steps = 50,
                        tol = 1e-10,
                        max_iter_inner = 200) {

  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))

  theta <- theta_init

  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )

  for (s in seq_len(n_steps)) {
    bw <- bw_grid[s]

    for (it in seq_len(max_iter_inner)) {
      pi_old     <- theta$pi
      lambda_old <- theta$lambda

      # E-step (r=1 고정)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = 1)

      # M-step (barrier on pi)
      theta <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )

      d_pi  <- abs(theta$pi - pi_old)
      d_lam <- abs(theta$lambda - lambda_old)

      trace <- rbind(trace, data.frame(
        step = s, inner = it, bw = bw,
        pi = theta$pi, lambda = theta$lambda,
        d_pi = d_pi, d_lam = d_lam
      ))

      if (max(d_pi, d_lam) < tol) break
    }
  }

  list(theta = theta, trace = trace, bw_grid = bw_grid)
}


# --- you provided (keep as-is) -------------------------------------------
zip_dkl_std_post <- function(x, theta0, theta1) {
  z0 <- (x == 0L)
  if (!any(z0)) return(0)

  g0 <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1 <- zip_estep_annealed(x[z0], pi = theta1$pi, lambda = theta1$lambda, r = 1)

  sum(g0 * (log(g0) - log(g1)) + (1 - g0) * (log(1 - g0) - log(1 - g1)))
}


zip_delta_dkl <- function(x, theta0, theta1, r) {
  z0 <- (x == 0L)
  if (!any(z0)) return(0)

  g_r0 <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = r)
  g0   <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1   <- zip_estep_annealed(x[z0], pi = theta1$pi, lambda = theta1$lambda, r = 1)

  sum(g_r0 * (log(g0) - log(g1)) + (1 - g_r0) * (log(1 - g0) - log(1 - g1)))
}

zip_delta_dkl_debug <- function(x, theta0, theta1, r, eps = 0, verbose = TRUE) {
  z0 <- (x == 0L)
  if (!any(z0)) return(list(value = 0, info = "no zeros"))

  x0 <- x[z0]

  g_r0 <- zip_estep_annealed(x0, pi = theta0$pi, lambda = theta0$lambda, r = r)
  g0   <- zip_estep_annealed(x0, pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1   <- zip_estep_annealed(x0, pi = theta1$pi, lambda = theta1$lambda, r = 1)

  # (선택) eps>0이면 클리핑해서 NaN을 없애고 값은 계속 계산 가능
  if (eps > 0) {
    clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)
    g_r0 <- clip01(g_r0); g0 <- clip01(g0); g1 <- clip01(g1)
  }

  # 로그 계산 전 점검
  bad <- list(
    g_r0_nonfinite = which(!is.finite(g_r0)),
    g0_nonfinite   = which(!is.finite(g0)),
    g1_nonfinite   = which(!is.finite(g1)),
    g0_le0         = which(g0 <= 0),
    g0_ge1         = which(g0 >= 1),
    g1_le0         = which(g1 <= 0),
    g1_ge1         = which(g1 >= 1),
    one_minus_g0_le0 = which((1 - g0) <= 0),
    one_minus_g1_le0 = which((1 - g1) <= 0)
  )

  # 실제로 log가 Inf/NaN을 만드는지 확인
  lg0  <- log(g0)
  lg1  <- log(g1)
  l1g0 <- log(1 - g0)
  l1g1 <- log(1 - g1)

  bad$log_g0_nonfinite   <- which(!is.finite(lg0))
  bad$log_g1_nonfinite   <- which(!is.finite(lg1))
  bad$log1mg0_nonfinite  <- which(!is.finite(l1g0))
  bad$log1mg1_nonfinite  <- which(!is.finite(l1g1))

  # 항별로도 NaN 생기는지 체크
  term1 <- g_r0 * (lg0 - lg1)
  term2 <- (1 - g_r0) * (l1g0 - l1g1)

  bad$term1_nonfinite <- which(!is.finite(term1))
  bad$term2_nonfinite <- which(!is.finite(term2))

  value <- sum(term1 + term2)

  if (verbose) {
    cat("---- zip_delta_dkl_debug ----\n")
    cat("n0 =", length(x0), "\n")
    cat("ranges:\n")
    cat("  g_r0:", sprintf("[%.3g, %.3g]", min(g_r0), max(g_r0)), "\n")
    cat("  g0  :", sprintf("[%.3g, %.3g]", min(g0), max(g0)), "\n")
    cat("  g1  :", sprintf("[%.3g, %.3g]", min(g1), max(g1)), "\n")
    cat("value =", value, "\n")

    # 문제가 있는 인덱스만 요약 출력
    show_idx <- unique(c(
      bad$g_r0_nonfinite, bad$g0_nonfinite, bad$g1_nonfinite,
      bad$g0_le0, bad$g0_ge1, bad$g1_le0, bad$g1_ge1,
      bad$log_g0_nonfinite, bad$log_g1_nonfinite,
      bad$log1mg0_nonfinite, bad$log1mg1_nonfinite,
      bad$term1_nonfinite, bad$term2_nonfinite
    ))
    show_idx <- show_idx[!is.na(show_idx)]

    if (length(show_idx) == 0) {
      cat("No non-finite issues detected.\n")
    } else {
      show_idx <- head(show_idx, 10)  # 너무 길면 10개만
      cat("Problem indices (first up to 10):", paste(show_idx, collapse=", "), "\n")
      df <- data.frame(
        i = show_idx,
        x0 = x0[show_idx],
        g_r0 = g_r0[show_idx],
        g0   = g0[show_idx],
        g1   = g1[show_idx],
        logg0 = lg0[show_idx],
        logg1 = lg1[show_idx],
        log1mg0 = l1g0[show_idx],
        log1mg1 = l1g1[show_idx],
        term1 = term1[show_idx],
        term2 = term2[show_idx]
      )
      print(df)
      cat("Bad sets sizes:\n")
      print(sapply(bad, length))
    }
    cat("-----------------------------\n")
  }

  list(value = value, bad = bad)
}


# DHEM with fixed schedules ------------------------------------------------
# - r_t : exp(seq(log(r_init), 0, length.out = n_steps))  (monotone up to 1)
# - bw_t: exp(seq(log(bw_init), log(bw_end), length.out = n_steps)) (monotone down if bw_end<bw_init)
# - At each step: E-step with annealed posterior (r_t),
#                 M-step with barrier on pi using bw_t

zip_dhem <- function(x,
                     theta_init,
                     r_init = 0.2,
                     bw_init = 5,
                     bw_end  = 0.05,
                     pi_min  = 0.005,
                     n_steps = 50,
                     tol = 1e-10,
                     max_iter_inner = 200) {

  # fixed schedules
  r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))                 # -> 1
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))      # bw_init -> bw_end

  theta <- theta_init

  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    r      = numeric(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )

  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]
    bw <- bw_grid[s]

    for (it in seq_len(max_iter_inner)) {
      pi_old     <- theta$pi
      lambda_old <- theta$lambda

      # E-step (annealed with r)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = r)

      # M-step (barrier on pi with bw)
      theta <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )

      d_pi  <- abs(theta$pi - pi_old)
      d_lam <- abs(theta$lambda - lambda_old)

      trace <- rbind(trace, data.frame(
        step = s, inner = it, r = r, bw = bw,
        pi = theta$pi, lambda = theta$lambda,
        d_pi = d_pi, d_lam = d_lam,
        stringsAsFactors = FALSE
      ))

      if (max(d_pi, d_lam) < tol) break
    }
  }

  list(theta = theta, trace = trace, r_grid = r_grid, bw_grid = bw_grid)
}


zip_adaptive_dhem <- function(x,
                              theta_init,
                              r_init  = 0.2,
                              bw_init = 1e-6,
                              pi_min  = 5e-6,
                              n_steps = 50,
                              eta = 0.1,
                              bw_rate = 0.1,
                              tol = 1e-10,
                              max_inner = 500,
                              verbose = TRUE
                            ) {

  r_grid   <- exp(seq(log(r_init), 0, length.out = n_steps))
  bw = bw_init
  theta <- theta_init

  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    r      = numeric(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  theta0 = list(pi=theta_init$pi,lambda=theta_init$lambda)

  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]

    theta_temp = theta0
    for (it in seq_len(max_inner)) {
      pi_old     <- theta_temp$pi
      lambda_old <- theta_temp$lambda
      # E-step (annealed with r)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = r)
      # M-step (barrier on pi with bw)
      theta1 <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )
      # print(theta1$pi)

      d_pi  <- abs(theta1$pi - theta_temp$pi )
      d_lam <- abs(theta1$lambda - theta_temp$lambda)

      if (max(d_pi, d_lam) < tol) {break
        }else{
        theta_temp=theta1
      }
    }

    dDKL   <- zip_delta_dkl(x, theta0 = theta0, theta1 = theta1, r = r)
    DKLstd <- zip_dkl_std_post(x, theta0 = theta0, theta1 = theta1)
    dB <- log(theta1$pi - pi_min) - log(theta0$pi - pi_min)

    Accept = FALSE
    if(is.na(dDKL-bw*dB)||dDKL-bw*dB<0){
      # Acceptance test start
      # Acc1
      if(dDKL<=eta*DKLstd){
        next
      }
      # Acc2
      if(eta*DKLstd<bw*abs(dB)){
        # bw = min(bw,eta*DKLstd/abs(dB))
        bw = min(bw,0.5*eta*DKLstd/abs(dB))
        Accept = FALSE
      }
    }else{
      Accept=TRUE
    }

    # parameter update
    if(Accept) theta0 = theta1
    #print(theta0$pi)
    trace <- rbind(trace, data.frame(
        step = s, inner = it, r = r, bw = bw,
        pi = theta0$pi, lambda = theta0$lambda,
        d_pi = d_pi, d_lam = d_lam,
        stringsAsFactors = FALSE
      ))

  }

  list(theta = theta0, bw = bw)
}


zip_diff <- function(x, pi, lambda, r = 1) {
  pi <- as.numeric(pi)[1]
  lambda <- as.numeric(lambda)[1]
  x <- as.numeric(x)

  gamma <- zip_estep_annealed(x, pi = pi, lambda = lambda, r = r)

  A <- sum(1 - gamma)
  B <- sum(gamma)

  A / pi - B / (1 - pi)
}

zip_bw_init_pi <- function(x, pi_init, lambda_init,
                           r_init = 1,
                           p_min = 0.005,
                           tau = 0.01) {
  pi_init <- as.numeric(pi_init)[1]
  lambda_init <- as.numeric(lambda_init)[1]
  x <- as.numeric(x)


  bw_raw <- tau * abs(zip_diff(x, pi_init, lambda_init, r = r_init)) *
    (pi_init - p_min)

  bw_raw
}
