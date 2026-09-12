## ================================================================
## DHEM.R
##
## Reproduction script for:
##   "Dual Homotopy Framework for Constrained EM Algorithm"
##
## This single script reproduces every figure and table in the paper's
## "Numerical Experiments" section:
##   1) Gaussian Mixture Model (GMM) simulation
##        -> Table: tab_gmm_err_mean_sd.csv          (Table "tab:gmm_err_mean_sd")
##   2) Zero-Inflated Poisson (ZIP) simulation
##        -> Figure: ZIPstandardEMresult.png          (fig:ZIP_stamdEMresult)
##        -> Figure: ZIPadpativeDHEMresult.png         (fig:ZIP_adaptiveDHEMresult)
##        -> Table:  tab_zip_err_mean_sd.csv           (Table "tab:zip_err_mean_sd")
##   3) Weibull Mixture Model (WMM) real-data analysis (bathtub curve)
##        -> Table:  tab_converged_params_math_header.csv (Table "tab:converged_params_math_header")
##        -> Figure: wmm_DAEM.png, wmm_BM.png, wmm_DHEM.png, wmm_adapDHEM.png (fig:trace-all)
##
## All model-fitting / simulation functions live in DHEM_functions.R
## (sourced below) and are used here UNMODIFIED from the original
## project code (gmm_functions.R, wmm_functions.R, ZIP_functions.R).
## Only this orchestration script was written to organize existing
## code into one reproducible entry point.
##
## Note on stochastic reproducibility: the original simulation code
## does not fix a global random seed (each replication draws its own
## seed via sample.int()), so re-running this script will not
## reproduce the paper's numbers bit-for-bit, but will reproduce them
## distributionally (same data-generating process, same sample sizes,
## same number of replications).
## ================================================================

## ---------------------------------------------------------------
## Simulation cache switch
## ---------------------------------------------------------------
## USE_SIM_CACHE = FALSE -> re-run the (randomized) simulations
##                          and overwrite the cached .rds files in DATA/
## USE_SIM_CACHE = TRUE  -> skip the simulations and load the previously
##                          saved results from DATA/ instead
USE_SIM_CACHE <- TRUE

## ---------------------------------------------------------------
## 0) Setup
## ---------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "clue", "patchwork",
                    "survival", "purrr", "scales", "parallel")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))

source("DHEM_functions.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(parallel)

fig_dir <- "figure"
tab_dir <- "table"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
if (!dir.exists(tab_dir)) dir.create(tab_dir)


data_dir <- "DATA"
if (!dir.exists(data_dir)) dir.create(data_dir)

gmm_cache_file <- file.path(data_dir, "gmm_sim_all.rds")
zip_cache_file <- file.path(data_dir, "zip_sim_cache.rds")

## method-name display order used throughout the paper's tables
gmm_zip_method_order <- c("AdaptiveDHEM", "Barrier", "DAEM", "DHEM", "EM")
method_display_name  <- c(AdaptiveDHEM = "Adaptive DHEM", Barrier = "Barrier",
                           DAEM = "DAEM", DHEM = "DHEM", EM = "EM")


## =================================================================
## 1) Gaussian Mixture Model (GMM) simulation
##    -> Table tab:gmm_err_mean_sd
## =================================================================

theta_true_gmm <- list(
  pi = c(0.2, 0.3, 0.5),
  mu = list(c(-1, 1, 2), c(1, 1, 0.5), c(2, 0, -2)),
  Sigma = list(
    matrix(c(
      1,     0.30, -0.20,
      0.30,  1,     0.10,
      -0.20, 0.10,  1
    ), 3, 3, byrow = TRUE),
    matrix(c(
      2,     0.18, -0.25,
      0.18,  0.10,  0.06,
      -0.25, 0.06,  0.50
    ), 3, 3, byrow = TRUE),
    matrix(c(
      0.50, -0.08,  0.22,
      -0.08,  0.10, -0.12,
      0.22, -0.12,  2
    ), 3, 3, byrow = TRUE)
  )
)

n_gmm       <- 100
K_gmm       <- 3
tol_limit   <- 1e-6
r_init_gmm  <- 0.1
n_steps_gmm <- 100
bw_init_gmm <- 1e-1
bw_end_gmm  <- 1e-10
delta_gmm   <- 1
max_iter_gmm <- 1000
eta_gmm     <- 0.5

## Number of Monte-Carlo replications. The paper reports M = 500; the
## original script used M = 100 for a faster run. Increase to 500 to
## match the paper exactly (runtime scales roughly linearly with M).
M_gmm <- 500

if (USE_SIM_CACHE) {
  gmm_sim_all <- readRDS(gmm_cache_file)
} else {
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("n_gmm", "K_gmm", "theta_true_gmm", "tol_limit", "r_init_gmm",
                       "n_steps_gmm", "bw_init_gmm", "bw_end_gmm", "delta_gmm",
                       "max_iter_gmm", "eta_gmm", "M_gmm"))
  invisible(clusterEvalQ(cl, { source("DHEM_functions.R") }))

  gmm_sim_list <- parLapply(cl, seq_len(M_gmm), function(m) {
    one_sim_summary(simulation_num = m, n = n_gmm, K = K_gmm, theta_true = theta_true_gmm,
                     tol_limit = tol_limit, r_init = r_init_gmm, n_steps = n_steps_gmm,
                     bw_init = bw_init_gmm, bw_end = bw_end_gmm, delta = delta_gmm,
                     max_iter = max_iter_gmm, eta = eta_gmm)
  })
  stopCluster(cl)

  gmm_sim_all <- do.call(rbind, gmm_sim_list)
  saveRDS(gmm_sim_all, gmm_cache_file)
}

gmm_err_table <- gmm_sim_all %>%
  filter(succ == 1) %>%
  group_by(method, para) %>%
  summarise(value = sprintf("%.3f (%.3f)", mean(err), sd(err)), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = para, values_from = value)

gmm_success_table <- gmm_sim_all %>%
  group_by(method) %>%
  summarise(success = sum(succ == 1, na.rm = TRUE) / n(), .groups = "drop")

tab_gmm_err_mean_sd <- gmm_err_table %>%
  left_join(gmm_success_table, by = "method") %>%
  mutate(method = factor(method, levels = gmm_zip_method_order)) %>%
  arrange(method) %>%
  mutate(method = method_display_name[as.character(method)]) %>%
  select(Method = method, pi, mu, Sigma, success)

print(tab_gmm_err_mean_sd)
write.csv(tab_gmm_err_mean_sd,
          file.path(tab_dir, "tab_gmm_err_mean_sd.csv"), row.names = FALSE)
writeLines(capture.output(print(tab_gmm_err_mean_sd)),
           file.path(tab_dir, "tab_gmm_err_mean_sd.txt"))


## =================================================================
## 2) Zero-Inflated Poisson (ZIP) simulation
##    -> Figure fig:ZIP_EM_vs_adaptiveDHEM (ZIPstandardEMresult.png,
##       ZIPadpativeDHEMresult.png), Table tab:zip_err_mean_sd
## =================================================================

n_zip           <- 1e4
pi_true_zip     <- 0.99
lambda_true_zip <- 0.4

pi_init_zip     <- 0.6
lambda_init_zip <- 1
theta_init_zip  <- list(pi = pi_init_zip, lambda = lambda_init_zip)
r_init_zip      <- 0.1
p_min_zip       <- 0.5
n_step_zip      <- 100
bw_end_zip      <- 1e-8

M_zip <- 100  # paper: M = 100 replications, n = 10,000

## ---- 2a) Table: averaged L1 error (pi, lambda) across methods ----
zip_sim_one_rep <- function(i) {
  x <- rzip(n_zip, pi = pi_true_zip, lambda = lambda_true_zip)

  bw_init_zip <- min(zip_bw_init_pi(x, pi_init_zip, lambda_init_zip, r_init_zip, p_min_zip, tau = 0.1))

  safe_fit <- function(expr) {
    out <- tryCatch(expr, error = function(e) NULL)
    if (is.null(out) || is.null(out$theta) || is.null(out$theta$pi) || is.null(out$theta$lambda)) {
      list(succ = 0L, theta = list(pi = NA_real_, lambda = NA_real_))
    } else {
      list(succ = 1L, theta = out$theta)
    }
  }

  res_em       <- safe_fit(zip_em(x, pi_init = pi_init_zip))
  res_daem     <- safe_fit(zip_daem(x, theta_init = theta_init_zip, r_init = r_init_zip, n_steps = n_step_zip))
  res_barrier  <- safe_fit(zip_barrier(x, theta_init = theta_init_zip, n_steps = n_step_zip,
                                        bw_init = bw_init_zip, bw_end = bw_end_zip, pi_min = p_min_zip))
  res_dhem     <- safe_fit(zip_dhem(x, theta_init = theta_init_zip, r_init = r_init_zip,
                                     bw_init = bw_init_zip, bw_end = bw_end_zip, pi_min = p_min_zip,
                                     n_steps = n_step_zip))
  res_adapDhem <- safe_fit(zip_adaptive_dhem(
    x = x, theta_init = theta_init_zip, r_init = r_init_zip, bw_init = bw_init_zip,
    pi_min = p_min_zip, n_steps = n_step_zip, bw_rate = 0.9, tol = 1e-10,
    max_inner = 1e4, verbose = FALSE
  ))

  rbind(
    data.frame(simulationNum = i, method = "EM", succ = res_em$succ,
               pi_err = ifelse(res_em$succ == 1L, abs_err(res_em$theta$pi, pi_true_zip), NA_real_),
               lambda_err = ifelse(res_em$succ == 1L, abs_err(res_em$theta$lambda, lambda_true_zip), NA_real_)),
    data.frame(simulationNum = i, method = "DAEM", succ = res_daem$succ,
               pi_err = ifelse(res_daem$succ == 1L, abs_err(res_daem$theta$pi, pi_true_zip), NA_real_),
               lambda_err = ifelse(res_daem$succ == 1L, abs_err(res_daem$theta$lambda, lambda_true_zip), NA_real_)),
    data.frame(simulationNum = i, method = "Barrier", succ = res_barrier$succ,
               pi_err = ifelse(res_barrier$succ == 1L, abs_err(res_barrier$theta$pi, pi_true_zip), NA_real_),
               lambda_err = ifelse(res_barrier$succ == 1L, abs_err(res_barrier$theta$lambda, lambda_true_zip), NA_real_)),
    data.frame(simulationNum = i, method = "DHEM", succ = res_dhem$succ,
               pi_err = ifelse(res_dhem$succ == 1L, abs_err(res_dhem$theta$pi, pi_true_zip), NA_real_),
               lambda_err = ifelse(res_dhem$succ == 1L, abs_err(res_dhem$theta$lambda, lambda_true_zip), NA_real_)),
    data.frame(simulationNum = i, method = "AdaptiveDHEM", succ = res_adapDhem$succ,
               pi_err = ifelse(res_adapDhem$succ == 1L, abs_err(res_adapDhem$theta$pi, pi_true_zip), NA_real_),
               lambda_err = ifelse(res_adapDhem$succ == 1L, abs_err(res_adapDhem$theta$lambda, lambda_true_zip), NA_real_))
  )
}

## Standard EM, M_zip replications
zip_em_one_rep <- function(m) {
  x <- rzip(n_zip, pi = pi_true_zip, lambda = lambda_true_zip)
  res_em <- zip_em(x, pi_init = pi_init_zip, lambda_init = lambda_init_zip)
  data.frame(sim = m, pi_hat = res_em$theta$pi, lambda_hat = res_em$theta$lambda)
}

## Adaptive DHEM, M_zip replications
zip_adap_one_rep <- function(m) {
  x <- rzip(n_zip, pi = pi_true_zip, lambda = lambda_true_zip)
  res <- zip_adaptive_dhem(
    x = x, theta_init = theta_init_zip, r_init = r_init_zip, bw_init = bw_init_zip,
    pi_min = p_min_zip, n_steps = n_step_zip, bw_rate = 0.9, tol = 1e-10,
    max_inner = 1e4, verbose = FALSE
  )
  data.frame(sim = m, pi_hat = res$theta$pi, lambda_hat = res$theta$lambda)
}

if (USE_SIM_CACHE) {
  zip_cache <- readRDS(zip_cache_file)
  sim_df  <- zip_cache$sim_df
  em_df   <- zip_cache$em_df
  adap_df <- zip_cache$adap_df
} else {
  ## bw_init for the Adaptive DHEM boxplot run below: the original script
  ## reused whatever bw_init_value happened to be left over from the last
  ## replication of its table loop. Here that loop runs in parallel, so no
  ## such leftover state exists; a representative bw_init is computed
  ## explicitly once from a fresh draw instead.
  x_bw_probe <- rzip(n_zip, pi = pi_true_zip, lambda = lambda_true_zip)
  bw_init_zip <- min(zip_bw_init_pi(x_bw_probe, pi_init_zip, lambda_init_zip, r_init_zip, p_min_zip, tau = 0.1))

  ## ZIP fitting is not parallelized in the original scripts (unlike GMM),
  ## and each replication re-fits 5 methods on n = 10,000 points, so a
  ## shared parallel cluster is used here for the table replications
  ## and the two boxplot replications below. This changes only *how*
  ## the (unmodified) functions are dispatched, not what they compute.
  n_cores_zip <- max(1, parallel::detectCores() - 1)
  cl_zip <- makeCluster(n_cores_zip)
  clusterExport(cl_zip, c("n_zip", "pi_true_zip", "lambda_true_zip",
                           "pi_init_zip", "lambda_init_zip", "theta_init_zip",
                           "r_init_zip", "p_min_zip", "n_step_zip", "bw_end_zip",
                           "bw_init_zip"))
  invisible(clusterEvalQ(cl_zip, { source("DHEM_functions.R") }))

  sim_df  <- do.call(rbind, parLapply(cl_zip, seq_len(M_zip), zip_sim_one_rep))
  em_df   <- bind_rows(parLapply(cl_zip, seq_len(M_zip), zip_em_one_rep))
  adap_df <- bind_rows(parLapply(cl_zip, seq_len(M_zip), zip_adap_one_rep)) %>%
    filter(is.finite(pi_hat), is.finite(lambda_hat))

  stopCluster(cl_zip)

  saveRDS(list(sim_df = sim_df, em_df = em_df, adap_df = adap_df), zip_cache_file)
}

## ---- 2a) Table: averaged L1 error (pi, lambda) across methods ----
tab_zip_err_mean_sd <- sim_df %>%
  group_by(method) %>%
  summarise(
    pi     = sprintf("%.3f (%.3f)", mean(pi_err, na.rm = TRUE), sd(pi_err, na.rm = TRUE)),
    lambda = sprintf("%.3f (%.3f)", mean(lambda_err, na.rm = TRUE), sd(lambda_err, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = gmm_zip_method_order)) %>%
  arrange(method) %>%
  mutate(Method = method_display_name[as.character(method)]) %>%
  select(Method, pi, lambda)

print(tab_zip_err_mean_sd)
write.csv(tab_zip_err_mean_sd,
          file.path(tab_dir, "tab_zip_err_mean_sd.csv"), row.names = FALSE)
writeLines(capture.output(print(tab_zip_err_mean_sd)),
           file.path(tab_dir, "tab_zip_err_mean_sd.txt"))

## ---- 2b) Figure: boxplots of pi/lambda estimates, standard EM vs adaptive DHEM ----

p_pi_em <- ggplot(em_df, aes(x = 1, y = pi_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = pi_true_zip, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "π (EM estimate)", x = NULL, y = "") +
  coord_cartesian(ylim = c(0.85, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

p_lam_em <- ggplot(em_df, aes(x = 1, y = lambda_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = lambda_true_zip, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "λ (EM estimate)", x = NULL, y = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

fig_zip_em <- (p_pi_em | p_lam_em)
ggsave(file.path(fig_dir, "ZIPstandardEMresult.png"), fig_zip_em, width = 8, height = 5, dpi = 150)

p_pi_adap <- ggplot(adap_df, aes(x = 1, y = pi_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = pi_true_zip, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "π (Adaptive DHEM estimate)", x = NULL, y = "") +
  coord_cartesian(ylim = c(0.85, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

p_lam_adap <- ggplot(adap_df, aes(x = 1, y = lambda_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = lambda_true_zip, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "λ (Adaptive DHEM estimate)", x = NULL, y = "") +
  coord_cartesian(ylim = c(0, 0.8)) +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

fig_zip_adap <- (p_pi_adap | p_lam_adap)
ggsave(file.path(fig_dir, "ZIPadpativeDHEMresult.png"), fig_zip_adap, width = 8, height = 5, dpi = 150)


## =================================================================
## 3) Weibull Mixture Model (WMM): real-data analysis (bathtub curve)
##    -> Table tab:converged_params_math_header, Figure fig:trace-all
## =================================================================

df_wmm <- read.table(file.path("DATA", "Aarest_data.txt"), header = TRUE)

maxGEMiter_wmm <- 1e6
nsteps_wmm     <- 100
r_init_wmm     <- 0.1
r_end_wmm      <- 1
bw_end_wmm     <- 1e-8
eta_wmm        <- 0.1
errtol_wmm     <- 1e-10

K_wmm <- 3
pi_init_wmm   <- rep(1 / K_wmm, K_wmm)
beta_init_wmm <- c(0.2, 1, 5)
lambda_init_wmm <- wmm_lambda_init(df_wmm$time, df_wmm$event, beta_init_wmm, ratio1 = 0.2, ratio3 = 0.8)
theta_init_wmm <- list(beta = beta_init_wmm, pi = pi_init_wmm, lambda = lambda_init_wmm)

bw_init_wmm <- wmm_bw_init_beta(
  df = df_wmm, pi_init = pi_init_wmm, lambda_init = lambda_init_wmm,
  beta_init = beta_init_wmm, r_init = r_init_wmm, tau = 0.1
)

fit_EM       <- wmm_EM(df_wmm, theta = theta_init_wmm, maxGEMiter = maxGEMiter_wmm, tol = errtol_wmm, verbose = FALSE)
fit_DAEM     <- wmm_DAEM(df_wmm, theta = theta_init_wmm, maxGEMiter = maxGEMiter_wmm, nsteps = nsteps_wmm,
                          r_init = r_init_wmm, r_end = r_end_wmm, tol = errtol_wmm, verbose = FALSE)
fit_BM       <- wmm_BM(df_wmm, theta = theta_init_wmm, maxGEMiter = maxGEMiter_wmm, nsteps = nsteps_wmm,
                        bw_init = bw_init_wmm, bw_end = bw_end_wmm, tol = errtol_wmm, verbose = FALSE)
fit_DHEM     <- wmm_DHEM(df_wmm, theta = theta_init_wmm, maxGEMiter = maxGEMiter_wmm, nsteps = nsteps_wmm,
                          r_init = r_init_wmm, r_end = r_end_wmm, bw_init = bw_init_wmm, bw_end = bw_end_wmm,
                          tol = errtol_wmm, verbose = FALSE)
fit_adapDHEM <- wmm_DHEM_adaptive(df_wmm, theta = theta_init_wmm, maxGEMiter = maxGEMiter_wmm, nsteps = nsteps_wmm,
                                   r_init = r_init_wmm, r_end = r_end_wmm, bw_init = bw_init_wmm, eta = eta_wmm,
                                   tol = errtol_wmm, verbose = FALSE)

df_EM       <- make_result_df(fit_EM, "EM")
df_DAEM     <- make_result_df(fit_DAEM, "DAEM")
df_BM       <- make_result_df(fit_BM, "Barrier")
df_DHEM     <- make_result_df(fit_DHEM, "DHEM")
df_adapDHEM <- make_result_df(fit_adapDHEM, "Adaptive DHEM")

## ---- Table: converged parameter estimates by method ----
## Use each method's LAST trace row (final convergence point at r=r_end,
## bw=bw_end) rather than best_row()'s minimum-|gradient| row: for DHEM the
## r-annealing schedule makes the GEM path pass through a transient
## near-zero-gradient region mid-trace (see wmm_DHEM.png), so best_row()
## can pick that transient point instead of where the method actually ends up.
last_row <- function(df) {
  df[nrow(df), c("method","pi1","pi2","pi3","beta1","beta3",
                 "lambda1","lambda2","lambda3","dQbeta1","dQbeta3")]
}

tab_converged_params_math_header <- do.call(rbind, lapply(
  list(df_EM, df_DAEM, df_BM, df_DHEM, df_adapDHEM),
  last_row
))

print(tab_converged_params_math_header)
write.csv(tab_converged_params_math_header,
          file.path(tab_dir, "tab_converged_params_math_header.csv"), row.names = FALSE)
writeLines(capture.output(print(tab_converged_params_math_header)),
           file.path(tab_dir, "tab_converged_params_math_header.txt"))

## ---- Figure: parameter and gradient traces by method ----
r_init <- r_init_wmm  # plot_beta_trace()/plot_dQbeta_trace() reference r_init from the calling scope

fig_wmm_DAEM <- (plot_beta_trace(df_DAEM, title = expression(DAEM * ": " * beta * " trace")) |
                    plot_dQbeta_trace(df_DAEM, title = expression(DAEM * ": " * nabla ~ beta * " trace"))) +
  plot_annotation(title = "DAEM") & theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(file.path(fig_dir, "wmm_DAEM.png"), fig_wmm_DAEM, width = 10, height = 8, dpi = 150)

fig_wmm_BM <- (plot_beta_trace_barrier(df_BM, title = expression(Barrier * " " * method * ": " * beta * " trace")) |
                 plot_dQbeta_trace_barrier(df_BM, title = expression(Barrier * " " * method * ": " * nabla ~ beta * " trace"))) +
  plot_annotation(title = "Barrier method") & theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(file.path(fig_dir, "wmm_BM.png"), fig_wmm_BM, width = 10, height = 8, dpi = 150)

fig_wmm_DHEM <- (plot_beta_trace(df_DHEM, title = expression(DHEM * ": " * beta * " trace")) |
                   plot_dQbeta_trace(df_DHEM, title = expression(DHEM * ": " * nabla ~ beta * " trace"))) +
  plot_annotation(title = "DHEM") & theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(file.path(fig_dir, "wmm_DHEM.png"), fig_wmm_DHEM, width = 10, height = 8, dpi = 150)

fig_wmm_adapDHEM <- (plot_beta_trace(df_adapDHEM, title = expression(adapDHEM * ": " * beta * " trace")) |
                        plot_dQbeta_trace(df_adapDHEM, title = expression(adapDHEM * ": " * nabla ~ beta * " trace"))) +
  plot_annotation(title = "Adaptive DHEM") & theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(file.path(fig_dir, "wmm_adapDHEM.png"), fig_wmm_adapDHEM, width = 10, height = 8, dpi = 150)


## ================================================================
## Done. Outputs written to ./figure and ./table
## ================================================================
message("Reproduction complete. Tables written to '", tab_dir, "/', figures written to '", fig_dir, "/'.")
