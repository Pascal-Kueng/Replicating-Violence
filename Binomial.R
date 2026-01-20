# ============================================================
# FULL ANALYSIS SCRIPT (Model 3 replication + diagnostics
# + alternative Binomial / Beta-binomial models + DHARMa)
#
# - Keeps your Stata-replicating intreg likelihood EXACTLY the same
# - Adds:
#   (A) Full diagnostics for the Tobit/intreg model
#   (B) Sensitivity sweep over censoring point cpoint
#   (C) Binomial GLM with cluster-robust SE
#   (D) Beta-binomial (glmmTMB) + corrected zero-prob simulation
#   (E) DHARMa diagnostics (if installed)
#   (F) Optional: random intercept / zero-inflation variants (try)
# ============================================================

# ----------------
# Packages
# ----------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# Optional packages (used if available)
have_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# ----------------
# Load data
# ----------------
df_raw <- read_excel("S1_file_combined.xlsx")

# ----------------
# Build analysis dataset (Model 3)
# ----------------
dat <- df_raw %>%
  transmute(
    tr_nop = as.numeric(tr_nop),   # trauma share in [0,1]
    ncases = as.numeric(ncases),
    
    co2LF  = as.factor(co2LF),     # clustering var
    co2    = as.factor(co2),       # FE dummies
    
    T0n_   = as.numeric(`T0n_`),
    T2be   = as.numeric(T2be),
    T3bm   = as.numeric(T3bm),
    T4bl   = as.numeric(T4bl),
    T5i_   = as.numeric(`T5i_`),
    
    battle = as.numeric(battle),
    royal  = as.numeric(royal),
    urban  = as.numeric(urban)
  ) %>%
  filter(
    !is.na(tr_nop), !is.na(ncases), !is.na(co2LF), !is.na(co2),
    !is.na(T0n_), !is.na(T2be), !is.na(T3bm), !is.na(T4bl), !is.na(T5i_),
    !is.na(battle), !is.na(royal), !is.na(urban),
    ncases > 0
  )

stopifnot(nrow(dat) == 82)

# Match Stata xi: base category if present (edit if your base differs)
if ("ir" %in% levels(dat$co2)) dat$co2 <- relevel(dat$co2, ref = "ir")

# ----------------
# Design matrix for Model 3
# ----------------
X <- model.matrix(
  ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2,
  data = dat
)
K <- ncol(X)
N <- nrow(dat)

# ----------------
# Tobit/intreg outcome construction (unchanged)
# ----------------
y <- ifelse(dat$tr_nop > 0, log(dat$tr_nop), NA_real_)
is_cens <- dat$tr_nop == 0

# Stata-normalized aweights: sum(a) = N  (unchanged)
a <- dat$ncases
a <- a * N / sum(a)

cluster <- dat$co2LF
stopifnot(length(cluster) == N)

# ============================================================
# PART 1: Tobit/intreg Model 3 (cpoint default = -3)
# ============================================================

cpoint <- -3

loglik <- function(par) {
  beta <- par[1:K]
  lns  <- par[K + 1]
  s <- exp(lns)
  xb <- as.vector(X %*% beta)
  
  ll <- numeric(N)
  
  # Uncensored
  idx_u <- which(!is_cens)
  r <- y[idx_u] - xb[idx_u]
  ll[idx_u] <- -0.5 * ( a[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a[idx_u]) )
  
  # Left-censored: y* <= cpoint
  idx_c <- which(is_cens)
  t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / s
  ll[idx_c] <- pnorm(t, log.p = TRUE)
  
  sum(ll)
}

start <- c(rep(0, K), 0)
fit <- optim(start, fn = function(p) -loglik(p), method = "BFGS", hessian = TRUE)
if (fit$convergence != 0) stop("optim did not converge")

theta <- fit$par
beta_hat <- theta[1:K]
lnsigma_hat <- theta[K + 1]
sigma_hat <- exp(lnsigma_hat)

xb <- as.vector(X %*% beta_hat)
sigma <- sigma_hat

# Scores U (unchanged)
U <- matrix(0, nrow = N, ncol = K + 1)
colnames(U) <- c(colnames(X), "lnsigma")

# Uncensored scores
idx_u <- which(!is_cens)
r <- y[idx_u] - xb[idx_u]
U[idx_u, 1:K]   <- (a[idx_u] * r / sigma^2) * X[idx_u, , drop = FALSE]
U[idx_u, K + 1] <- -1 + (a[idx_u] * r^2 / sigma^2)

# Censored scores
idx_c <- which(is_cens)
t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / sigma
log_phi <- dnorm(t, log = TRUE)
log_Phi <- pnorm(t, log.p = TRUE)
lambda <- exp(pmin(log_phi - log_Phi, 700))

U[idx_c, 1:K]   <- -(lambda * sqrt(a[idx_c]) / sigma) * X[idx_c, , drop = FALSE]
U[idx_c, K + 1] <- -(lambda * t)

# Cluster-robust sandwich (unchanged)
bread <- solve(fit$hessian)
Sg <- rowsum(U, cluster)
meat <- crossprod(Sg)
V_CR0 <- bread %*% meat %*% bread
G <- nrow(Sg)
V_stata <- (G/(G - 1)) * V_CR0

Vb <- V_stata[1:K, 1:K, drop = FALSE]
se <- sqrt(diag(Vb))
zv <- beta_hat / se
pv <- 2 * pnorm(-abs(zv))

out_stata <- data.frame(
  Estimate     = beta_hat,
  `Std. Error` = se,
  `z value`    = zv,
  `Pr(>|z|)`   = pv,
  row.names    = colnames(X)
)

cat("\n====================\nModel 3 (Tobit/intreg) coefficients\n====================\n")
print(out_stata)
cat("\nSigma:\n")
print(c(sigma = sigma_hat, lnsigma = lnsigma_hat))

# ============================================================
# PART 2: Tobit/intreg diagnostics (baseline cpoint)
# ============================================================

# Fit stats
ll_hat <- loglik(theta)
k_par <- K + 1
AIC_t <- -2 * ll_hat + 2 * k_par
BIC_t <- -2 * ll_hat + log(N) * k_par

# Wald chi2 of slopes (robust VCE)
idx_slope <- which(colnames(X) != "(Intercept)")
b_slope <- beta_hat[idx_slope]
V_slope <- Vb[idx_slope, idx_slope, drop = FALSE]
Wald_chi2 <- as.numeric(t(b_slope) %*% solve(V_slope) %*% b_slope)
df_wald <- length(idx_slope)
p_wald <- pchisq(Wald_chi2, df = df_wald, lower.tail = FALSE)

# Censoring probabilities (model implied)
sd_i <- sigma_hat / sqrt(a)
p_cens_hat <- pnorm((cpoint - xb) / sd_i)

obs_cens_unw  <- mean(is_cens)
pred_cens_unw <- mean(p_cens_hat)
obs_cens_w    <- weighted.mean(is_cens, w = a)
pred_cens_w   <- weighted.mean(p_cens_hat, w = a)

# Weighted Brier score for censoring indicator (0/1)
brier_w <- weighted.mean((as.numeric(is_cens) - p_cens_hat)^2, w = a)

# RMSE among uncensored (descriptive)
res_u <- y[!is_cens] - xb[!is_cens]
rmse_u <- sqrt(weighted.mean(res_u^2, w = a[!is_cens]))

# Null model + McFadden pseudo-R2 (recomputed for this cpoint)
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
if (fit0$convergence != 0) stop("Null (intercept-only) optim did not converge")
ll_null <- loglik0(fit0$par)
pseudoR2_mcfadden <- 1 - (ll_hat / ll_null)

# Surprise ranking for censoring indicator
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

fit_stats <- data.frame(
  N = N,
  uncensored = sum(!is_cens),
  left_censored = sum(is_cens),
  clusters = G,
  logLik = ll_hat,
  AIC = AIC_t,
  BIC = BIC_t,
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

cat("\n====================\nTobit/intreg diagnostics (baseline cpoint)\n====================\n")
print(fit_stats)

cat("\nTop weighted surprises (baseline cpoint):\n")
print(top_surprise)

# Plots (optional)
par(mfrow = c(1, 2))
plot(xb[!is_cens], res_u,
     xlab = "Fitted xb (latent mean)",
     ylab = "Residual (y - xb) [uncensored]",
     main = "Uncensored residuals vs fitted")
abline(h = 0, lty = 2)

boxplot(p_cens_hat ~ is_cens,
        names = c("Uncensored", "Censored"),
        ylab = "Predicted P(censored)",
        main = "Predicted censoring probability")
par(mfrow = c(1, 1))

# ============================================================
# PART 3: Sensitivity sweep over cpoint (same full pipeline core)
# ============================================================

run_intreg_for_cpoint <- function(cpoint, start_par) {
  
  loglik_c <- function(par) {
    beta <- par[1:K]
    lns  <- par[K + 1]
    s <- exp(lns)
    xb <- as.vector(X %*% beta)
    
    ll <- numeric(N)
    
    idx_u <- which(!is_cens)
    r <- y[idx_u] - xb[idx_u]
    ll[idx_u] <- -0.5 * ( a[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a[idx_u]) )
    
    idx_c <- which(is_cens)
    t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / s
    ll[idx_c] <- pnorm(t, log.p = TRUE)
    
    sum(ll)
  }
  
  fit <- optim(start_par, fn = function(p) -loglik_c(p), method = "BFGS", hessian = TRUE)
  if (fit$convergence != 0) stop("optim did not converge at cpoint=", cpoint)
  
  theta <- fit$par
  beta_hat <- theta[1:K]
  sigma_hat <- exp(theta[K + 1])
  xb_hat <- as.vector(X %*% beta_hat)
  
  ll_hat <- loglik_c(theta)
  AIC <- -2 * ll_hat + 2 * (K + 1)
  BIC <- -2 * ll_hat + log(N) * (K + 1)
  
  sd_i <- sigma_hat / sqrt(a)
  p_cens_hat <- pnorm((cpoint - xb_hat) / sd_i)
  
  brier_w <- weighted.mean((as.numeric(is_cens) - p_cens_hat)^2, w = a)
  
  p_unc_hat <- 1 - p_cens_hat
  surprise <- ifelse(is_cens, -log(p_cens_hat), -log(p_unc_hat))
  w_surprise <- a * surprise
  
  data.frame(
    cpoint = cpoint,
    logLik = ll_hat,
    AIC = AIC,
    BIC = BIC,
    obs_cens_unw = mean(is_cens),
    pred_cens_unw = mean(p_cens_hat),
    obs_cens_w = weighted.mean(is_cens, w = a),
    pred_cens_w = weighted.mean(p_cens_hat, w = a),
    brier_w = brier_w,
    sigma = sigma_hat,
    max_w_surprise = max(w_surprise)
  )
}

c_grid <- c(-5, -4.5, -4, -3.5, -3, -2.5)
compare <- do.call(rbind, {
  start_par <- theta  # warm-start from baseline solution
  out <- vector("list", length(c_grid))
  for (j in seq_along(c_grid)) {
    out[[j]] <- run_intreg_for_cpoint(c_grid[j], start_par = start_par)
    # update warm start by re-optimizing at this cpoint (best practice)
    # (re-using last theta is good enough; we keep it simple)
  }
  out
})

cat("\n====================\nCensor-point sensitivity (Tobit/intreg)\n====================\n")
print(compare[order(compare$BIC), ])

# ============================================================
# PART 4: Binomial / Beta-binomial models using ncases as trials
# ============================================================

# ---- Reconstruct integer counts k from proportion p and trials n ----
dat2 <- dat %>%
  mutate(
    n = as.integer(round(ncases)),
    p = tr_nop
  )

stopifnot(all(dat2$n > 0))
stopifnot(all(dat2$p >= 0 & dat2$p <= 1))

# Deterministic rounding (simple)
dat2 <- dat2 %>%
  mutate(
    k_raw = p * n,
    k = as.integer(round(k_raw)),
    k = pmin(pmax(k, 0L), n)
  )

recon_check <- dat2 %>%
  mutate(diff = k_raw - k) %>%
  summarise(
    max_abs_diff = max(abs(diff)),
    mean_abs_diff = mean(abs(diff)),
    share_exact_integer = mean(abs(diff) < 1e-8)
  )

cat("\n====================\nCount reconstruction check (k ≈ p*n)\n====================\n")
print(recon_check)

# ---- Binomial GLM (logit) + cluster-robust SEs ----
f_bin <- cbind(k, n - k) ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2

m_bin <- glm(f_bin, data = dat2, family = binomial(link = "logit"))

cat("\n====================\nBinomial GLM (logit)\n====================\n")
print(summary(m_bin)$coef)

if (have_pkg("clubSandwich") && have_pkg("lmtest")) {
  library(clubSandwich)
  library(lmtest)
  
  V_CR <- vcovCR(m_bin, cluster = dat2$co2LF, type = "CR0")
  cat("\nCluster-robust (CR0) inference for binomial GLM:\n")
  print(coeftest(m_bin, vcov. = V_CR))
} else {
  cat("\nNOTE: Install clubSandwich + lmtest for cluster-robust SEs in GLM.\n")
}

# Overdispersion check (binomial)
overdisp_ratio <- sum(residuals(m_bin, type = "pearson")^2) / df.residual(m_bin)
cat("\nOverdispersion ratio (Pearson / df): ", overdisp_ratio, "\n", sep = "")

# ---- Beta-binomial via glmmTMB ----

  library(glmmTMB)
  
  m_bb <- glmmTMB(f_bin, data = dat2, family = betabinomial(link = "logit"))
  cat("\n====================\nBeta-binomial (glmmTMB)\n====================\n")
  print(summary(m_bb))
  
  cat("\nAIC/BIC comparison (binomial vs beta-binomial):\n")
  print(AIC(m_bin, m_bb))
  print(BIC(m_bin, m_bb))
  
  # ---- Zero probability calibration (binomial, closed form) ----
  p_hat_bin <- predict(m_bin, type = "response")
  p0_hat_bin <- (1 - p_hat_bin) ^ dat2$n
  
  cat("\nZero-prob calibration (binomial):\n")
  print(c(
    obs_zero_unw = mean(dat2$k == 0),
    pred_zero_unw = mean(p0_hat_bin),
    obs_zero_w_n = weighted.mean(dat2$k == 0, w = dat2$n),
    pred_zero_w_n = weighted.mean(p0_hat_bin, w = dat2$n)
  ))
  
  # ---- Zero probability calibration (beta-binomial, simulation) ----
  set.seed(1)
  sim_obj <- simulate(m_bb, nsim = 500)
  
  # IMPORTANT: simulate(m_bb) returns a data.frame with 500 columns;
  # each column is an N x 2 matrix: [k, n-k]. We ONLY want k (first column).
  sim_k_mat <- do.call(cbind, lapply(sim_obj, function(m) m[, 1]))
  stopifnot(nrow(sim_k_mat) == nrow(dat2), ncol(sim_k_mat) == 500)
  
  p0_hat_bb <- rowMeans(sim_k_mat == 0)
  stopifnot(length(p0_hat_bb) == nrow(dat2))
  
  cat("\nZero-prob calibration (beta-binomial; simulated):\n")
  print(c(
    obs_zero_unw = mean(dat2$k == 0),
    pred_zero_unw = mean(p0_hat_bb),
    obs_zero_w_n = weighted.mean(dat2$k == 0, w = dat2$n),
    pred_zero_w_n = weighted.mean(p0_hat_bb, w = dat2$n)
  ))
  
  # ---- DHARMa diagnostics (if installed) ----
  if (have_pkg("DHARMa")) {
    library(DHARMa)
    
    cat("\n====================\nDHARMa diagnostics: beta-binomial\n====================\n")
    res_bb <- simulateResiduals(m_bb, plot = TRUE)
    print(testDispersion(res_bb))
    print(testZeroInflation(res_bb))
  } else {
    cat("\nNOTE: Install DHARMa for residual/dispersion/zero-inflation diagnostics.\n")
  }
  
  # ---- Optional sensitivity models (try/catch so script won't break) ----
  cat("\n====================\nOptional models (may or may not improve)\n====================\n")
  
  # Random intercept for co2LF (alternative to cluster-robust)
  m_bb_re <- try(
    glmmTMB(update(f_bin, . ~ . + (1 | co2LF)), data = dat2, family = betabinomial(link = "logit")),
    silent = TRUE
  )
  if (!inherits(m_bb_re, "try-error")) {
    cat("\nBeta-binomial + random intercept (1|co2LF):\n")
    print(AIC(m_bb, m_bb_re))
    print(BIC(m_bb, m_bb_re))
  } else {
    cat("\nCould not fit beta-binomial + (1|co2LF) (convergence or unsupported).\n")
  }
  
  # Zero-inflated binomial (often works; ZI beta-binomial may or may not)
  m_bin_zi <- try(
    glmmTMB(f_bin, ziformula = ~ 1, data = dat2, family = binomial(link = "logit")),
    silent = TRUE
  )
  if (!inherits(m_bin_zi, "try-error")) {
    cat("\nZero-inflated binomial (ziformula ~1):\n")
    print(AIC(m_bin, m_bin_zi))
    print(BIC(m_bin, m_bin_zi))
    if (have_pkg("DHARMa")) {
      cat("\nDHARMa diagnostics: ZI binomial\n")
      res_zi <- simulateResiduals(m_bin_zi, plot = TRUE)
      print(testDispersion(res_zi))
      print(testZeroInflation(res_zi))
    }
  } else {
    cat("\nCould not fit zero-inflated binomial.\n")
  }
  
  m_bb_zi <- try(
    glmmTMB(f_bin, ziformula = ~ 1, data = dat2, family = betabinomial(link = "logit")),
    silent = TRUE
  )
  if (!inherits(m_bb_zi, "try-error")) {
    cat("\nZero-inflated beta-binomial (ziformula ~1):\n")
    print(AIC(m_bb, m_bb_zi))
    print(BIC(m_bb, m_bb_zi))
    if (have_pkg("DHARMa")) {
      cat("\nDHARMa diagnostics: ZI beta-binomial\n")
      res_bbzi <- simulateResiduals(m_bb_zi, plot = TRUE)
      print(testDispersion(res_bbzi))
      print(testZeroInflation(res_bbzi))
    }
  } else {
    cat("\nCould not fit zero-inflated beta-binomial (common; depends on build/data).\n")
  }# ============================================================
  # ADD-ON DIAGNOSTICS: cluster-robust inference for glmmTMB
  # + predictive checks (overall + by period)
  # ============================================================
  
  # --- helper: safe package attach ---
  have_pkg <- function(p) requireNamespace(p, quietly = TRUE)
  
  # ------------------------------------------------------------
  # 1) Cluster-robust SEs for glmmTMB beta-binomial (and variants)
  # ------------------------------------------------------------
  if (exists("m_bb") && have_pkg("clubSandwich")) {
    library(clubSandwich)
    
    cat("\n====================\nCluster-robust inference: beta-binomial (glmmTMB)\n====================\n")
    
    # CR2 is often preferred for small samples; CR0 is the classic sandwich.
    # We'll try CR2 first, then fall back to CR0 if needed.
    V_bb <- try(vcovCR(m_bb, cluster = dat2$co2LF, type = "CR2"), silent = TRUE)
    if (inherits(V_bb, "try-error")) {
      cat("vcovCR CR2 failed; trying CR0...\n")
      V_bb <- vcovCR(m_bb, cluster = dat2$co2LF, type = "CR0")
    }
    
    # coef_test gives robust tests; "naive-t" uses a t reference (often better w/ few clusters)
    # If that errors, fall back to z.
    ct_bb <- try(coef_test(m_bb, vcov = V_bb, test = "naive-t"), silent = TRUE)
    if (inherits(ct_bb, "try-error")) {
      ct_bb <- coef_test(m_bb, vcov = V_bb, test = "naive-z")
    }
    print(ct_bb)
    
    # If you also fit these, you get robust inference too:
    for (nm in c("m_bb_re", "m_bb_zi", "m_bin_zi")) {
      if (exists(nm)) {
        mod <- get(nm)
        cat("\n--------------------\nCluster-robust inference for:", nm, "\n--------------------\n")
        V_ <- try(vcovCR(mod, cluster = dat2$co2LF, type = "CR2"), silent = TRUE)
        if (inherits(V_, "try-error")) V_ <- vcovCR(mod, cluster = dat2$co2LF, type = "CR0")
        
        ct_ <- try(coef_test(mod, vcov = V_, test = "naive-t"), silent = TRUE)
        if (inherits(ct_, "try-error")) ct_ <- coef_test(mod, vcov = V_, test = "naive-z")
        print(ct_)
      }
    }
  } else {
    cat("\nNOTE: clubSandwich not available or m_bb not in environment; skipping glmmTMB cluster-robust inference.\n")
  }
  
  # ------------------------------------------------------------
  # 2) Define "period" labels from the time dummies
  #    (assumes mutually exclusive indicators; otherwise it's still usable)
  # ------------------------------------------------------------
  dat2 <- dat2 %>%
    mutate(
      period = case_when(
        T0n_ == 1 ~ "T0n_",
        T2be == 1 ~ "T2be",
        T3bm == 1 ~ "T3bm",
        T4bl == 1 ~ "T4bl",
        T5i_ == 1 ~ "T5i_",
        TRUE ~ "BASE"
      )
    )
  
  # ------------------------------------------------------------
  # 3) Common-scale predictions: predicted mean share and predicted zero rate
  # ------------------------------------------------------------
  
  # --- Observed ---
  obs_share <- dat2$p
  obs_zero  <- (dat2$k == 0)
  
  # --- Binomial GLM predictions (mean share and P(k=0)) ---
  if (exists("m_bin")) {
    p_hat_bin <- predict(m_bin, type = "response")     # E[k/n] under binomial
    p0_hat_bin <- (1 - p_hat_bin) ^ dat2$n             # P(k=0 | n,p)
  } else {
    p_hat_bin <- rep(NA_real_, N)
    p0_hat_bin <- rep(NA_real_, N)
  }
  
  # --- Beta-binomial predictions (mean share and P(k=0) via simulation) ---
  # Mean share:
  if (exists("m_bb")) {
    p_hat_bb <- as.numeric(predict(m_bb, type = "response"))  # E[k/n] (conditional mean)
    
    # P(k=0) via simulation (CORRECT: use successes column only)
    set.seed(1)
    sim_obj <- simulate(m_bb, nsim = 500)
    sim_k_mat <- do.call(cbind, lapply(sim_obj, function(m) m[, 1]))  # successes only
    stopifnot(nrow(sim_k_mat) == nrow(dat2), ncol(sim_k_mat) == 500)
    p0_hat_bb <- rowMeans(sim_k_mat == 0)
  } else {
    p_hat_bb <- rep(NA_real_, N)
    p0_hat_bb <- rep(NA_real_, N)
  }
  
  # --- Tobit/intreg implied predictions on original SHARE scale ---
  # We'll do a simulation-based predictive mean share and P(zero) on the same "0 vs positive" definition.
  # Here: predict zero if latent y* <= cpoint; share = 0 if censored, else exp(y*).
  if (exists("xb") && exists("sigma_hat") && exists("a") && exists("cpoint")) {
    set.seed(1)
    sd_i <- sigma_hat / sqrt(a)
    S <- 5000
    # simulate latent y* for each obs: y* ~ N(xb, sd_i)
    # use matrix operations: N x S
    ystar <- matrix(rnorm(N * S), nrow = N, ncol = S)
    ystar <- xb + sd_i * ystar
    
    share_sim <- ifelse(ystar <= cpoint, 0, exp(ystar))
    share_hat_tobit <- rowMeans(share_sim)
    p0_hat_tobit <- rowMeans(ystar <= cpoint)
  } else {
    share_hat_tobit <- rep(NA_real_, N)
    p0_hat_tobit <- rep(NA_real_, N)
  }
  
  # ------------------------------------------------------------
  # 4) Summary tables: overall + by period, weighted/unweighted
  #    Use weights = n (trials) as the natural weighting for proportions.
  #    Also show weights = a (your Stata-normalized analytic weights) for comparability.
  # ------------------------------------------------------------
  
  summ_overall <- data.frame(
    metric = c("mean_share_unw", "mean_share_w_n", "mean_share_w_a",
               "zero_rate_unw", "zero_rate_w_n", "zero_rate_w_a"),
    observed = c(
      mean(obs_share),
      weighted.mean(obs_share, w = dat2$n),
      weighted.mean(obs_share, w = a),
      mean(obs_zero),
      weighted.mean(obs_zero, w = dat2$n),
      weighted.mean(obs_zero, w = a)
    ),
    tobit = c(
      mean(share_hat_tobit),
      weighted.mean(share_hat_tobit, w = dat2$n),
      weighted.mean(share_hat_tobit, w = a),
      mean(p0_hat_tobit),
      weighted.mean(p0_hat_tobit, w = dat2$n),
      weighted.mean(p0_hat_tobit, w = a)
    ),
    binomial = c(
      mean(p_hat_bin),
      weighted.mean(p_hat_bin, w = dat2$n),
      weighted.mean(p_hat_bin, w = a),
      mean(p0_hat_bin),
      weighted.mean(p0_hat_bin, w = dat2$n),
      weighted.mean(p0_hat_bin, w = a)
    ),
    betabinomial = c(
      mean(p_hat_bb),
      weighted.mean(p_hat_bb, w = dat2$n),
      weighted.mean(p_hat_bb, w = a),
      mean(p0_hat_bb),
      weighted.mean(p0_hat_bb, w = dat2$n),
      weighted.mean(p0_hat_bb, w = a)
    )
  )
  
  cat("\n====================\nPredictive checks (overall)\n====================\n")
  print(summ_overall)
  
  summ_by_period <- dat2 %>%
    transmute(
      period,
      n = n,
      w_a = a,
      obs_share = obs_share,
      obs_zero = obs_zero,
      tobit_share = share_hat_tobit,
      tobit_p0 = p0_hat_tobit,
      bin_share = p_hat_bin,
      bin_p0 = p0_hat_bin,
      bb_share = p_hat_bb,
      bb_p0 = p0_hat_bb
    ) %>%
    group_by(period) %>%
    summarise(
      n_obs = n(),
      # observed
      obs_share_unw = mean(obs_share),
      obs_share_w_n = weighted.mean(obs_share, w = n),
      obs_share_w_a = weighted.mean(obs_share, w = w_a),
      obs_zero_unw  = mean(obs_zero),
      obs_zero_w_n  = weighted.mean(obs_zero, w = n),
      obs_zero_w_a  = weighted.mean(obs_zero, w = w_a),
      # tobit
      tobit_share_unw = mean(tobit_share),
      tobit_share_w_n = weighted.mean(tobit_share, w = n),
      tobit_share_w_a = weighted.mean(tobit_share, w = w_a),
      tobit_p0_unw    = mean(tobit_p0),
      tobit_p0_w_n    = weighted.mean(tobit_p0, w = n),
      tobit_p0_w_a    = weighted.mean(tobit_p0, w = w_a),
      # binomial
      bin_share_unw = mean(bin_share),
      bin_share_w_n = weighted.mean(bin_share, w = n),
      bin_share_w_a = weighted.mean(bin_share, w = w_a),
      bin_p0_unw    = mean(bin_p0),
      bin_p0_w_n    = weighted.mean(bin_p0, w = n),
      bin_p0_w_a    = weighted.mean(bin_p0, w = w_a),
      # beta-binomial
      bb_share_unw = mean(bb_share),
      bb_share_w_n = weighted.mean(bb_share, w = n),
      bb_share_w_a = weighted.mean(bb_share, w = w_a),
      bb_p0_unw    = mean(bb_p0),
      bb_p0_w_n    = weighted.mean(bb_p0, w = n),
      bb_p0_w_a    = weighted.mean(bb_p0, w = w_a)
    )
  
  cat("\n====================\nPredictive checks by period\n====================\n")
  print(summ_by_period)
  
  # ------------------------------------------------------------
  # 5) Zero calibration tables (deciles) for each model
  # ------------------------------------------------------------
  make_calibration <- function(p0_hat, obs_zero, w, name) {
    cuts <- unique(quantile(p0_hat, probs = seq(0, 1, 0.1), na.rm = TRUE))
    bin <- cut(p0_hat, breaks = cuts, include.lowest = TRUE)
    out <- data.frame(bin = bin, p0_hat = p0_hat, obs_zero = obs_zero, w = w) %>%
      group_by(bin) %>%
      summarise(
        n = n(),
        pred_unw = mean(p0_hat),
        obs_unw  = mean(obs_zero),
        pred_w   = weighted.mean(p0_hat, w),
        obs_w    = weighted.mean(obs_zero, w)
      )
    out$model <- name
    out
  }
  
  cal_all <- bind_rows(
    make_calibration(p0_hat_tobit, obs_zero, w = dat2$n, name = "tobit"),
    make_calibration(p0_hat_bin,   obs_zero, w = dat2$n, name = "binomial"),
    make_calibration(p0_hat_bb,    obs_zero, w = dat2$n, name = "betabinomial")
  )
  
  cat("\n====================\nZero calibration by decile (weights = n)\n====================\n")
  print(cal_all)
  
  # ------------------------------------------------------------
  # 6) Optional quick plots: observed vs predicted by period (mean share + zero rate)
  # ------------------------------------------------------------
  op <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  
  # Mean share (weighted by n)
  plot(
    summ_by_period$obs_share_w_n, summ_by_period$bb_share_w_n,
    xlab = "Observed mean share (weighted by n)",
    ylab = "Predicted mean share (beta-binomial, weighted by n)",
    main = "Mean share by period (beta-binomial)",
    pch = 19
  )
  abline(0, 1, lty = 2)
  text(summ_by_period$obs_share_w_n, summ_by_period$bb_share_w_n, labels = summ_by_period$period, pos = 4, cex = 0.8)
  
  plot(
    summ_by_period$obs_share_w_n, summ_by_period$tobit_share_w_n,
    xlab = "Observed mean share (weighted by n)",
    ylab = "Predicted mean share (Tobit sim, weighted by n)",
    main = "Mean share by period (Tobit)",
    pch = 19
  )
  abline(0, 1, lty = 2)
  text(summ_by_period$obs_share_w_n, summ_by_period$tobit_share_w_n, labels = summ_by_period$period, pos = 4, cex = 0.8)
  
  # Zero rate (weighted by n)
  plot(
    summ_by_period$obs_zero_w_n, summ_by_period$bb_p0_w_n,
    xlab = "Observed zero rate (weighted by n)",
    ylab = "Predicted P(zero) (beta-binomial, weighted by n)",
    main = "Zero rate by period (beta-binomial)",
    pch = 19
  )
  abline(0, 1, lty = 2)
  text(summ_by_period$obs_zero_w_n, summ_by_period$bb_p0_w_n, labels = summ_by_period$period, pos = 4, cex = 0.8)
  
  plot(
    summ_by_period$obs_zero_w_n, summ_by_period$tobit_p0_w_n,
    xlab = "Observed zero rate (weighted by n)",
    ylab = "Predicted P(zero) (Tobit, weighted by n)",
    main = "Zero rate by period (Tobit)",
    pch = 19
  )
  abline(0, 1, lty = 2)
  text(summ_by_period$obs_zero_w_n, summ_by_period$tobit_p0_w_n, labels = summ_by_period$period, pos = 4, cex = 0.8)
  
  par(op)
  
  cat("\nADD-ON CHECKS COMPLETE.\n")
  