
# ============================================================
# Full pipeline re-run for a grid of censoring points
# (ONLY cpoint changes; everything else identical)
# ============================================================

source("Model_3.R")

run_full_pipeline_for_cpoint <- function(cpoint, start_par = NULL, seed_rq = 1) {
  
  # ---- loglik for this cpoint (identical structure) ----
  loglik_c <- function(par) {
    beta <- par[1:K]
    lns  <- par[K + 1]
    s <- exp(lns)
    xb <- as.vector(X %*% beta)
    
    ll <- numeric(N)
    
    # Uncensored
    idx_u <- which(!is_cens)
    r <- y[idx_u] - xb[idx_u]
    ll[idx_u] <- -0.5 * ( a[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a[idx_u]) )
    
    # Left-censored
    idx_c <- which(is_cens)
    t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / s
    ll[idx_c] <- pnorm(t, log.p = TRUE)
    
    sum(ll)
  }
  
  # ---- MLE ----
  if (is.null(start_par)) start_par <- c(rep(0, K), 0)
  fit <- optim(start_par, fn = function(p) -loglik_c(p), method = "BFGS", hessian = TRUE)
  if (fit$convergence != 0) stop("optim did not converge at cpoint = ", cpoint)
  
  theta <- fit$par
  beta_hat <- theta[1:K]
  lnsigma_hat <- theta[K + 1]
  sigma_hat <- exp(lnsigma_hat)
  xb <- as.vector(X %*% beta_hat)
  
  # ---- Scores U (same as baseline, but this cpoint) ----
  U <- matrix(0, nrow = N, ncol = K + 1)
  colnames(U) <- c(colnames(X), "lnsigma")
  
  # Uncensored scores
  idx_u <- which(!is_cens)
  r <- y[idx_u] - xb[idx_u]
  U[idx_u, 1:K]   <- (a[idx_u] * r / sigma_hat^2) * X[idx_u, , drop = FALSE]
  U[idx_u, K + 1] <- -1 + (a[idx_u] * r^2 / sigma_hat^2)
  
  # Censored scores
  idx_c <- which(is_cens)
  t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / sigma_hat
  log_phi <- dnorm(t, log = TRUE)
  log_Phi <- pnorm(t, log.p = TRUE)
  lambda <- exp(pmin(log_phi - log_Phi, 700))
  
  U[idx_c, 1:K]   <- -(lambda * sqrt(a[idx_c]) / sigma_hat) * X[idx_c, , drop = FALSE]
  U[idx_c, K + 1] <- -(lambda * t)
  
  # ---- Cluster-robust sandwich (same) ----
  bread <- solve(fit$hessian)
  Sg <- rowsum(U, cluster)
  meat <- crossprod(Sg)
  V_CR0 <- bread %*% meat %*% bread
  G <- nrow(Sg)
  V_stata <- (G/(G - 1)) * V_CR0
  
  Vb <- V_stata[1:K, 1:K, drop = FALSE]
  se <- sqrt(diag(Vb))
  z  <- beta_hat / se
  p  <- 2 * pnorm(-abs(z))
  
  out_stata <- data.frame(
    Estimate     = beta_hat,
    `Std. Error` = se,
    `z value`    = z,
    `Pr(>|z|)`   = p,
    row.names    = colnames(X)
  )
  
  # ---- Fit stats + censoring probs ----
  ll_hat <- loglik_c(theta)
  k_par <- K + 1
  AIC <- -2 * ll_hat + 2 * k_par
  BIC <- -2 * ll_hat + log(N) * k_par
  
  sd_i <- sigma_hat / sqrt(a)
  p_cens_hat <- pnorm((cpoint - xb) / sd_i)
  
  obs_cens_unw  <- mean(is_cens)
  pred_cens_unw <- mean(p_cens_hat)
  obs_cens_w    <- weighted.mean(is_cens, w = a)
  pred_cens_w   <- weighted.mean(p_cens_hat, w = a)
  
  brier_w <- weighted.mean((as.numeric(is_cens) - p_cens_hat)^2, w = a)
  
  # RMSE among uncensored (descriptive)
  res_u <- y[!is_cens] - xb[!is_cens]
  rmse_u <- sqrt(weighted.mean(res_u^2, w = a[!is_cens]))
  
  # ---- Null model + pseudo-R2 (must be recomputed per cpoint) ----
  X0 <- matrix(1, nrow = N, ncol = 1)
  colnames(X0) <- "(Intercept)"
  
  loglik0 <- function(par0) {
    b0 <- par0[1]
    lns0 <- par0[2]
    s0 <- exp(lns0)
    xb0 <- rep(b0, N)
    
    ll0 <- numeric(N)
    
    idx_u0 <- which(!is_cens)
    r0 <- y[idx_u0] - xb0[idx_u0]
    ll0[idx_u0] <- -0.5 * ( a[idx_u0] * (r0^2) / (s0^2) + log(2*pi) + 2*lns0 - log(a[idx_u0]) )
    
    idx_c0 <- which(is_cens)
    t0 <- (cpoint - xb0[idx_c0]) * sqrt(a[idx_c0]) / s0
    ll0[idx_c0] <- pnorm(t0, log.p = TRUE)
    
    sum(ll0)
  }
  
  fit0 <- optim(c(0, 0), fn = function(p) -loglik0(p), method = "BFGS", hessian = TRUE)
  if (fit0$convergence != 0) stop("Null optim did not converge at cpoint = ", cpoint)
  
  ll_null <- loglik0(fit0$par)
  pseudoR2_mcfadden <- 1 - (ll_hat / ll_null)
  
  # ---- Expected observed outcome E[y_obs|X] where y_obs = max(cpoint, y*) ----
  zc <- (xb - cpoint) / sd_i
  Ey_obs <- pnorm(zc) * xb + sd_i * dnorm(zc) + (1 - pnorm(zc)) * cpoint
  
  pred_summary <- data.frame(
    xb = xb,
    p_cens_hat = p_cens_hat,
    Ey_obs = Ey_obs
  )
  
  # ---- Calibration table (deciles) ----
  cuts <- quantile(p_cens_hat, probs = seq(0, 1, 0.1), na.rm = TRUE)
  cuts <- unique(cuts)
  bin <- cut(p_cens_hat, breaks = cuts, include.lowest = TRUE)
  
  cal <- data.frame(
    bin = bin,
    p_hat = p_cens_hat,
    cens = is_cens,
    w = a
  ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      n = dplyr::n(),
      pred_unw = mean(p_hat),
      obs_unw  = mean(cens),
      pred_w   = weighted.mean(p_hat, w),
      obs_w    = weighted.mean(cens, w)
    )
  
  # ---- Weighted "surprise" ranking (your narrative driver) ----
  p_unc_hat <- 1 - p_cens_hat
  surprise <- ifelse(is_cens, -log(p_cens_hat), -log(p_unc_hat))
  w_surprise <- a * surprise
  
  o_s <- order(w_surprise, decreasing = TRUE)
  top_surprise <- data.frame(
    row = o_s[1:min(10, N)],
    a = a[o_s[1:min(10, N)]],
    xb = xb[o_s[1:min(10, N)]],
    is_cens = is_cens[o_s[1:min(10, N)]],
    p_cens = p_cens_hat[o_s[1:min(10, N)]],
    surprise = surprise[o_s[1:min(10, N)]],
    w_surprise = w_surprise[o_s[1:min(10, N)]]
  )
  
  # ---- Outliers among uncensored (standardized residual) ----
  std_res_u <- res_u / sd_i[!is_cens]
  o_u <- order(abs(std_res_u), decreasing = TRUE)
  top_outliers <- data.frame(
    row = which(!is_cens)[o_u[1:min(10, length(o_u))]],
    xb = xb[!is_cens][o_u[1:min(10, length(o_u))]],
    y  = y[!is_cens][o_u[1:min(10, length(o_u))]],
    resid = res_u[o_u[1:min(10, length(o_u))]],
    std_resid = std_res_u[o_u[1:min(10, length(o_u))]],
    weight_a = a[!is_cens][o_u[1:min(10, length(o_u))]],
    Region_Period_Cluster_Unit = cluster[which(!is_cens)[o_u[1:min(10, length(o_u))]]]
  )
  
  # ---- Wald chi2 of slopes (exclude intercept) ----
  idx_slope <- which(colnames(X) != "(Intercept)")
  b_slope <- beta_hat[idx_slope]
  V_slope <- Vb[idx_slope, idx_slope, drop = FALSE]
  Wald_chi2 <- as.numeric(t(b_slope) %*% solve(V_slope) %*% b_slope)
  df_wald <- length(idx_slope)
  p_wald <- pchisq(Wald_chi2, df = df_wald, lower.tail = FALSE)
  
  fit_stats <- data.frame(
    cpoint = cpoint,
    N = N,
    uncensored = sum(!is_cens),
    left_censored = sum(is_cens),
    clusters = G,
    logLik = ll_hat,
    AIC = AIC,
    BIC = BIC,
    Wald_chi2 = Wald_chi2,
    df = df_wald,
    p_value = p_wald,
    obs_cens_unw = obs_cens_unw,
    pred_cens_unw = pred_cens_unw,
    obs_cens_w = obs_cens_w,
    pred_cens_w = pred_cens_w,
    brier_w = brier_w,
    rmse_uncensored = rmse_u,
    sigma = sigma_hat,
    ll_null = ll_null,
    pseudoR2_mcfadden = pseudoR2_mcfadden,
    max_w_surprise = max(w_surprise)
  )
  
  list(
    cpoint = cpoint,
    theta = theta,
    out = out_stata,
    fit_stats = fit_stats,
    pred_summary = pred_summary,
    cal = cal,
    top_surprise = top_surprise,
    top_outliers = top_outliers
  )
}

# ---- Choose a grid of censor points to try ----
c_grid <- c(-5, -4.5, -4, -3.5, -3, -2.5)

# ---- Run all models (warm-start each from previous theta for stability) ----
results <- vector("list", length(c_grid))

start_par <- c(beta_hat, log(sigma_hat))  # warm start from your baseline solution if available
for (j in seq_along(c_grid)) {
  results[[j]] <- run_full_pipeline_for_cpoint(c_grid[j], start_par = start_par)
  start_par <- results[[j]]$theta
}

# ---- Compare fits ----
compare <- do.call(rbind, lapply(results, function(r) r$fit_stats))
print(compare[order(compare$BIC), ])  # sort by BIC (or AIC)

# ---- If you want to inspect one cpoint in detail ----
# results[[which(c_grid == -4)]]$out
# results[[which(c_grid == -4)]]$cal
# results[[which(c_grid == -4)]]$top_surprise
# results[[which(c_grid == -4)]]$top_outliers


# ============================================================
# AGGREGATE RESULTS (FIXED)
# ============================================================

# 1. Fit Statistics Table (Sorted by cpoint)
fit_grid <- do.call(rbind, lapply(results, function(r) r$fit_stats))
fit_grid <- fit_grid[order(fit_grid$cpoint), ] 
print(fit_grid)

# 2. Coefficients Table (formatted with stars)
format_est_sig <- function(model_res) {
  # Extract Estimate (Column 1) and P-value (Column 4) by index to be safe
  est  <- model_res$out[, 1]  
  se   <- model_res$out[, 2]
  pval <- model_res$out[, 4]  
  
  # Define stars
  stars <- ifelse(pval < 0.01, "**",
                  ifelse(pval < 0.05, "*", ""))
  
  # Format string
  sprintf("%.3f%s (%.3f)", est, stars, se)
}

# Combine results into a matrix
est_matrix_sig <- do.call(cbind, lapply(results, format_est_sig))

# Set variable names
var_names <- rownames(results[[1]]$out)

# Set column names based on cpoint
colnames(est_matrix_sig) <- paste0("c=", sapply(results, function(r) r$cpoint))

# Combine into data frame
coef_grid_sig <- data.frame(
  Variable = var_names, 
  est_matrix_sig, 
  check.names = FALSE
)

# Load library (install if needed: install.packages("knitr"))
library(knitr)

# ============================================================
# OPTION 1: Static "Kable" Tables (Good for PDF/Word/Console)
# ============================================================

# 1. Coefficients Table
#    (Already formatted as strings with stars, so no rounding needed here)
print(kable(coef_grid_sig, 
            caption = "Table 1: Coefficients across censoring points (* p<0.05, ** p<0.01)",
            align = 'c'))

# 2. Fit Statistics Table
#    (Selecting key columns to avoid width issues, and rounding numbers)
key_cols <- c("cpoint", "logLik", "AIC", "BIC", "rmse_uncensored", 
              "pseudoR2_mcfadden", "obs_cens_w", "pred_cens_w")

print(kable(fit_grid[, key_cols], 
            digits = 3, 
            caption = "Table 2: Model Fit Statistics by Censoring Point",
            align = 'c'))


library(ggplot2)

# ============================================================
# 3. COEFFICIENT PLOT (Visual Sensitivity Analysis) - FIXED
# ============================================================

# 1. Extract data from the results list into a single long data frame
plot_data <- do.call(rbind, lapply(results, function(r) {
  df <- r$out
  df$Variable <- rownames(df)
  df$cpoint   <- r$cpoint
  
  # FIX: Rename the standard error column safely (handles "Std. Error", "Std..Error", etc.)
  # We look for the column index that contains "Std"
  se_col_idx <- grep("Std", colnames(df))
  colnames(df)[se_col_idx] <- "SE"
  
  return(df)
}))

# 2. Create the Confidence Intervals (95%)
# Now plot_data$SE exists, so this math will work
plot_data$lower <- plot_data$Estimate - 1.96 * plot_data$SE
plot_data$upper <- plot_data$Estimate + 1.96 * plot_data$SE

# 3. Plot
# We facet by Variable so each coefficient gets its own panel
p <- ggplot(plot_data, aes(x = cpoint, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_line(color = "blue") +
  geom_point(size = 1.5) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Sensitivity of Coefficients to Censoring Point",
    subtitle = "Shaded area = 95% Confidence Interval",
    x = "Censoring Point (c)",
    y = "Estimate (beta)"
  )

print(p)


report::report_packages()
