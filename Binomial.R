# ============================================================
# FULL ANALYSIS SCRIPT (Model 3 replication + diagnostics
# + alternative Binomial / Beta-binomial models + DHARMa
# + common-scale predictive checks + glmmTMB cluster-robust SEs)
#
# Key goals:
#  1) Tobit/intreg code is unchanged estimation-wise (Stata replica),
#     only Model 3 covariates + i.co2 added to X.
#  2) Add diagnostics + sensitivity on censoring point.
#  3) Add binomial + beta-binomial models on reconstructed counts.
#  4) Add DHARMa checks.
#  5) Add common-scale predictive checks and manual cluster-robust
#     inference for glmmTMB (conditional fixed effects only).
# ============================================================

# ----------------
# Packages
# ----------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

have_pkg <- function(p) requireNamespace(p, quietly = TRUE)

add_dharma_title <- function(model_label) {
  mtext(paste0("DHARMa residual diagnostics: ", model_label),
        side = 3, line = 1, cex = 1.0, font = 2)
}


# ----------------
# Load data
# ----------------
data_url <- "https://zenodo.org/records/8010025/files/S1_file_combined.xlsx?download=1"
tmp_file <- tempfile(fileext = ".xlsx")
download.file(data_url, destfile = tmp_file, mode = "wb", quiet = TRUE)
df_raw <- read_excel(tmp_file)
unlink(tmp_file)

# ----------------
# Build analysis dataset (Model 3)
# ----------------
dat <- df_raw %>%
  transmute(
    tr_nop = as.numeric(tr_nop),   # trauma share in [0,1]
    ncases = as.numeric(ncases),
    
    co2LF  = as.factor(co2LF),     # clustering var (as in your Model 1)
    co2    = as.factor(co2),       # i.co2 dummies
    
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

# Match Stata xi: base category for i.co2 (if relevant)
if ("ir" %in% levels(dat$co2)) dat$co2 <- relevel(dat$co2, ref = "ir")

# ----------------
# Design matrix for Model 3
# ----------------
X <- model.matrix(
  ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2,
  data = dat
)  # includes intercept

K <- ncol(X)
N <- nrow(dat)

# ----------------
# Tobit/intreg outcome construction (unchanged)
# ----------------
y <- ifelse(dat$tr_nop > 0, log(dat$tr_nop), NA_real_)
is_cens <- dat$tr_nop == 0

# Stata-normalized aweights: sum(a) = N (unchanged)
a <- dat$ncases
a <- a * N / sum(a)

cluster <- dat$co2LF
stopifnot(length(cluster) == N)

# ============================================================
# PART 1: Tobit/intreg Model 3 (baseline cpoint = -3)
# ============================================================

cpoint <- -3

# --- Log-likelihood: Stata intreg with aweights (sigma/sqrt(a)) ---
loglik <- function(par) {
  beta <- par[1:K]
  lns  <- par[K + 1]
  s <- exp(lns)
  xb <- as.vector(X %*% beta)
  
  ll <- numeric(N)
  
  # Uncensored: y observed exactly
  idx_u <- which(!is_cens)
  r <- y[idx_u] - xb[idx_u]
  ll[idx_u] <- -0.5 * ( a[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a[idx_u]) )
  
  # Left-censored: y* <= cpoint
  idx_c <- which(is_cens)
  t <- (cpoint - xb[idx_c]) * sqrt(a[idx_c]) / s
  ll[idx_c] <- pnorm(t, log.p = TRUE)
  
  sum(ll)
}

# --- Fit MLE ---
start <- c(rep(0, K), 0)  # betas=0, lnsigma=0
fit <- optim(start, fn = function(p) -loglik(p), method = "BFGS", hessian = TRUE)
if (fit$convergence != 0) stop("optim did not converge")

theta <- fit$par
beta_hat <- theta[1:K]
lnsigma_hat <- theta[K + 1]
sigma_hat <- exp(lnsigma_hat)

xb <- as.vector(X %*% beta_hat)
sigma <- sigma_hat

# --- Per-observation scores U (N x (K+1)) for the log-likelihood ---
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

# lambda = phi(t)/Phi(t), computed stably
log_phi <- dnorm(t, log = TRUE)
log_Phi <- pnorm(t, log.p = TRUE)
lambda <- exp(pmin(log_phi - log_Phi, 700))

U[idx_c, 1:K]   <- -(lambda * sqrt(a[idx_c]) / sigma) * X[idx_c, , drop = FALSE]
U[idx_c, K + 1] <- -(lambda * t)

# --- Cluster-robust sandwich variance ---
bread <- solve(fit$hessian)              # inverse Hessian of NEGATIVE log-likelihood
Sg <- rowsum(U, cluster)                # cluster-summed scores
meat <- crossprod(Sg)
V_CR0 <- bread %*% meat %*% bread

# --- Stata-style cluster small-sample correction for ML: multiply by G/(G-1) only ---
G <- nrow(Sg)
V_stata <- (G/(G - 1)) * V_CR0

# --- Output table (betas only) ---
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

# Model-implied censoring probs
sd_i <- sigma_hat / sqrt(a)
p_cens_hat <- pnorm((cpoint - xb) / sd_i)

obs_cens_unw  <- mean(is_cens)
pred_cens_unw <- mean(p_cens_hat)
obs_cens_w    <- weighted.mean(is_cens, w = a)
pred_cens_w   <- weighted.mean(p_cens_hat, w = a)

# Weighted Brier score for censoring indicator
brier_w <- weighted.mean((as.numeric(is_cens) - p_cens_hat)^2, w = a)

# RMSE among uncensored (descriptive)
res_u <- y[!is_cens] - xb[!is_cens]
rmse_u <- sqrt(weighted.mean(res_u^2, w = a[!is_cens]))

# Null model + McFadden pseudo-R2
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

# Optional diagnostic plots
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
# PART 3: Sensitivity sweep over censor point cpoint
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
  
  fit_ <- optim(start_par, fn = function(p) -loglik_c(p), method = "BFGS", hessian = TRUE)
  if (fit_$convergence != 0) stop("optim did not converge at cpoint=", cpoint)
  
  th <- fit_$par
  beta_h <- th[1:K]
  sigma_h <- exp(th[K + 1])
  xb_h <- as.vector(X %*% beta_h)
  
  ll_h <- loglik_c(th)
  AIC <- -2 * ll_h + 2 * (K + 1)
  BIC <- -2 * ll_h + log(N) * (K + 1)
  
  sd_i <- sigma_h / sqrt(a)
  p_cens <- pnorm((cpoint - xb_h) / sd_i)
  
  brier_w <- weighted.mean((as.numeric(is_cens) - p_cens)^2, w = a)
  
  p_unc <- 1 - p_cens
  surprise <- ifelse(is_cens, -log(p_cens), -log(p_unc))
  w_surprise <- a * surprise
  
  data.frame(
    cpoint = cpoint,
    logLik = ll_h,
    AIC = AIC,
    BIC = BIC,
    obs_cens_unw = mean(is_cens),
    pred_cens_unw = mean(p_cens),
    obs_cens_w = weighted.mean(is_cens, w = a),
    pred_cens_w = weighted.mean(p_cens, w = a),
    brier_w = brier_w,
    sigma = sigma_h,
    max_w_surprise = max(w_surprise)
  )
}

c_grid <- c(-5, -4.5, -4, -3.5, -3, -2.5)
compare <- do.call(rbind, {
  start_par <- theta
  out <- vector("list", length(c_grid))
  for (j in seq_along(c_grid)) {
    out[[j]] <- run_intreg_for_cpoint(c_grid[j], start_par = start_par)
  }
  out
})

cat("\n====================\nCensor-point sensitivity (Tobit/intreg)\n====================\n")
print(compare[order(compare$BIC), ])

# ============================================================
# PART 4: Binomial + Beta-binomial models on reconstructed counts
# ============================================================

dat2 <- dat %>%
  mutate(
    n = as.integer(round(ncases)),
    p = tr_nop
  )

stopifnot(all(dat2$n > 0))
stopifnot(all(dat2$p >= 0 & dat2$p <= 1))

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

# Binomial GLM (logit)
f_bin <- cbind(k, n - k) ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2
m_bin <- glm(f_bin, data = dat2, family = binomial(link = "logit"))

cat("\n====================\nBinomial GLM (logit): model-based SE\n====================\n")
print(summary(m_bin)$coef)

# Cluster-robust SE for binomial GLM
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

# Beta-binomial via glmmTMB + DHARMa + optional models
m_bb <- NULL
m_bb_re <- NULL
m_bb_zi <- NULL
m_bin_zi <- NULL

if (!have_pkg("glmmTMB")) {
  cat("\nNOTE: Install glmmTMB to fit beta-binomial / ZI models.\n")
} else {
  library(glmmTMB)
  
  m_bb <- glmmTMB(f_bin, data = dat2, family = betabinomial(link = "logit"))
  cat("\n====================\nBeta-binomial (glmmTMB)\n====================\n")
  print(summary(m_bb))
  
  cat("\nAIC/BIC comparison (binomial vs beta-binomial):\n")
  print(AIC(m_bin, m_bb))
  print(BIC(m_bin, m_bb))
  
  # Zero probability calibration (binomial, closed form)
  p_hat_bin <- predict(m_bin, type = "response")
  p0_hat_bin <- (1 - p_hat_bin) ^ dat2$n
  
  cat("\nZero-prob calibration (binomial):\n")
  print(c(
    obs_zero_unw   = mean(dat2$k == 0),
    pred_zero_unw  = mean(p0_hat_bin),
    obs_zero_w_n   = weighted.mean(dat2$k == 0, w = dat2$n),
    pred_zero_w_n  = weighted.mean(p0_hat_bin, w = dat2$n),
    obs_zero_w_a   = weighted.mean(dat2$k == 0, w = a),
    pred_zero_w_a  = weighted.mean(p0_hat_bin, w = a)
  ))
  
  # Zero probability calibration (beta-binomial, simulation; SUCCESS COUNTS ONLY)
  set.seed(1)
  sim_obj <- simulate(m_bb, nsim = 500)
  sim_k_mat <- do.call(cbind, lapply(sim_obj, function(m) m[, 1]))
  stopifnot(nrow(sim_k_mat) == nrow(dat2), ncol(sim_k_mat) == 500)
  p0_hat_bb <- rowMeans(sim_k_mat == 0)
  
  cat("\nZero-prob calibration (beta-binomial; simulated):\n")
  print(c(
    obs_zero_unw   = mean(dat2$k == 0),
    pred_zero_unw  = mean(p0_hat_bb),
    obs_zero_w_n   = weighted.mean(dat2$k == 0, w = dat2$n),
    pred_zero_w_n  = weighted.mean(p0_hat_bb, w = dat2$n),
    obs_zero_w_a   = weighted.mean(dat2$k == 0, w = a),
    pred_zero_w_a  = weighted.mean(p0_hat_bb, w = a)
  ))
  
  # DHARMa diagnostics
  if (have_pkg("DHARMa")) {
    library(DHARMa)
    cat("\n====================\nDHARMa diagnostics: beta-binomial\n====================\n")
    res_bb <- simulateResiduals(m_bb, plot = FALSE)
    plot(res_bb)
    add_dharma_title("Beta-binomial (glmmTMB)")
    print(testDispersion(res_bb))
    print(testZeroInflation(res_bb))
  } else {
    cat("\nNOTE: Install DHARMa for residual/dispersion/zero-inflation diagnostics.\n")
  }
  
  # Optional models (try; script will not break on failure)
  cat("\n====================\nOptional models (may or may not improve)\n====================\n")
  
  m_bb_re <- try(
    glmmTMB(update(f_bin, . ~ . + (1 | co2LF)), data = dat2, family = betabinomial(link = "logit")),
    silent = TRUE
  )
  if (!inherits(m_bb_re, "try-error")) {
    cat("\nBeta-binomial + random intercept (1|co2LF):\n")
    print(AIC(m_bb, m_bb_re))
    print(BIC(m_bb, m_bb_re))
  } else {
    m_bb_re <- NULL
    cat("\nCould not fit beta-binomial + (1|co2LF).\n")
  }
  
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
      res_zi <- simulateResiduals(m_bin_zi, plot = FALSE)
      plot(res_zi)
      add_dharma_title("Zero-inflated binomial (glmmTMB)")
      print(testDispersion(res_zi))
      print(testZeroInflation(res_zi))
    }
  } else {
    m_bin_zi <- NULL
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
      res_bbzi <- simulateResiduals(m_bb_zi, plot = FALSE)
      plot(res_bbzi)
      add_dharma_title("Zero-inflated beta-binomial (glmmTMB)")
      print(testDispersion(res_bbzi))
      print(testZeroInflation(res_bbzi))
    }
  } else {
    m_bb_zi <- NULL
    cat("\nCould not fit zero-inflated beta-binomial.\n")
  }
}

# ============================================================
# PART 5: ADD-ON CHECKS
#   A) Manual cluster-robust inference for glmmTMB (conditional FE)
#   B) Common-scale predictive checks (overall + by period)
#   C) Zero calibration by deciles
#   D) Optional plots
# ============================================================

# --- A) Manual cluster-robust inference for glmmTMB conditional FE ---
if (!is.null(m_bb) && have_pkg("clubSandwich")) {
  library(clubSandwich)
  
  robust_glmmTMB_cond <- function(mod, cluster_vec, type = c("CR2", "CR0")) {
    type <- match.arg(type)
    V_full <- try(vcovCR(mod, cluster = cluster_vec, type = type), silent = TRUE)
    if (inherits(V_full, "try-error")) return(NULL)
    
    b <- fixef(mod)$cond
    cn <- colnames(V_full)
    
    idx <- match(names(b), cn)
    if (anyNA(idx)) {
      idx2 <- match(paste0("cond_", names(b)), cn)
      if (!anyNA(idx2)) idx <- idx2
    }
    if (anyNA(idx)) return(NULL)
    
    V <- V_full[idx, idx, drop = FALSE]
    se <- sqrt(diag(V))
    z  <- as.numeric(b) / se
    p  <- 2 * pnorm(-abs(z))
    
    data.frame(
      Estimate = as.numeric(b),
      `CR SE`  = se,
      z = z,
      p = p,
      row.names = names(b)
    )
  }
  
  cat("\n====================\nBeta-binomial (glmmTMB) cluster-robust inference (conditional FE)\n====================\n")
  tab_bb <- robust_glmmTMB_cond(m_bb, dat2$co2LF, type = "CR2")
  if (is.null(tab_bb)) {
    cat("CR2 failed; trying CR0...\n")
    tab_bb <- robust_glmmTMB_cond(m_bb, dat2$co2LF, type = "CR0")
  }
  if (is.null(tab_bb)) {
    cat("Could not compute robust table for m_bb (name matching issue). Try: head(colnames(vcovCR(m_bb,...)))\n")
  } else {
    print(tab_bb)
  }
  
  for (nm in c("m_bb_re", "m_bb_zi", "m_bin_zi")) {
    mod <- get(nm)
    if (!is.null(mod)) {
      cat("\n--------------------\nCluster-robust (conditional FE) for:", nm, "\n--------------------\n")
      tab <- robust_glmmTMB_cond(mod, dat2$co2LF, type = "CR2")
      if (is.null(tab)) tab <- robust_glmmTMB_cond(mod, dat2$co2LF, type = "CR0")
      if (is.null(tab)) {
        cat("Could not compute robust table for ", nm, "\n", sep = "")
      } else {
        print(tab)
      }
    }
  }
} else {
  cat("\nNOTE: Skipping glmmTMB cluster-robust inference (need m_bb + clubSandwich).\n")
}

# --- B) Period labels for predictive checks ---
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

# --- Common-scale observed outcomes ---
obs_share <- dat2$p
obs_zero  <- (dat2$k == 0)

# --- Binomial predictions ---
p_hat_bin <- predict(m_bin, type = "response")
p0_hat_bin <- (1 - p_hat_bin) ^ dat2$n

# --- Beta-binomial predictions: mean share + P(zero) via simulation ---
p_hat_bb <- rep(NA_real_, N)
p0_hat_bb <- rep(NA_real_, N)
if (!is.null(m_bb)) {
  p_hat_bb <- as.numeric(predict(m_bb, type = "response"))
  
  set.seed(1)
  sim_obj <- simulate(m_bb, nsim = 500)
  sim_k_mat <- do.call(cbind, lapply(sim_obj, function(m) m[, 1]))
  stopifnot(nrow(sim_k_mat) == nrow(dat2), ncol(sim_k_mat) == 500)
  p0_hat_bb <- rowMeans(sim_k_mat == 0)
}

# --- Tobit/intreg predictions on original SHARE scale via simulation ---
share_hat_tobit <- rep(NA_real_, N)
p0_hat_tobit <- rep(NA_real_, N)
{
  set.seed(1)
  sd_i <- sigma_hat / sqrt(a)
  S <- 5000
  ystar <- matrix(rnorm(N * S), nrow = N, ncol = S)
  ystar <- xb + sd_i * ystar
  
  share_sim <- ifelse(ystar <= cpoint, 0, exp(ystar))
  share_hat_tobit <- rowMeans(share_sim)
  p0_hat_tobit <- rowMeans(ystar <= cpoint)
}

# --- Overall predictive checks (unweighted, weighted by n, weighted by a) ---
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
    mean(p_hat_bb, na.rm = TRUE),
    weighted.mean(p_hat_bb, w = dat2$n, na.rm = TRUE),
    weighted.mean(p_hat_bb, w = a, na.rm = TRUE),
    mean(p0_hat_bb, na.rm = TRUE),
    weighted.mean(p0_hat_bb, w = dat2$n, na.rm = TRUE),
    weighted.mean(p0_hat_bb, w = a, na.rm = TRUE)
  )
)

cat("\n====================\nPredictive checks (overall)\n====================\n")
print(summ_overall)

# --- Predictive checks by period (weighted by n) ---
summ_by_period <- dat2 %>%
  transmute(
    period,
    n = n,
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
    obs_share_w_n = weighted.mean(obs_share, w = n),
    tobit_share_w_n = weighted.mean(tobit_share, w = n),
    bin_share_w_n = weighted.mean(bin_share, w = n),
    bb_share_w_n = weighted.mean(bb_share, w = n, na.rm = TRUE),
    
    obs_zero_w_n = weighted.mean(obs_zero, w = n),
    tobit_p0_w_n = weighted.mean(tobit_p0, w = n),
    bin_p0_w_n = weighted.mean(bin_p0, w = n),
    bb_p0_w_n = weighted.mean(bb_p0, w = n, na.rm = TRUE)
  )

cat("\n====================\nPredictive checks by period (weighted by n)\n====================\n")
print(summ_by_period)

# --- C) Zero calibration by deciles (weights = n) ---
make_calibration <- function(p0_hat, obs_zero, w, name) {
  cuts <- unique(quantile(p0_hat, probs = seq(0, 1, 0.1), na.rm = TRUE))
  bin <- cut(p0_hat, breaks = cuts, include.lowest = TRUE)
  out <- data.frame(bin = bin, p0_hat = p0_hat, obs_zero = obs_zero, w = w) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      pred_w = weighted.mean(p0_hat, w, na.rm = TRUE),
      obs_w  = weighted.mean(obs_zero, w, na.rm = TRUE)
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

# --- D) Optional plots: observed vs predicted zero rate by period (with titles + subtitle) ---
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

add_subtitle <- function(txt) {
  mtext(txt, side = 3, line = 0.2, cex = 0.85)
}

good_line_txt <- "Good fit: points close to 45° line; no systematic over/under-prediction across periods."

# 1) Beta-binomial: zero rate by period
plot(
  summ_by_period$obs_zero_w_n, summ_by_period$bb_p0_w_n,
  xlab = "Observed zero rate (weighted by n)",
  ylab = "Predicted P(zero) (weighted by n)",
  main = "Beta-binomial (glmmTMB): Zero-rate calibration by period",
  pch = 19
)
add_subtitle(good_line_txt)
abline(0, 1, lty = 2)
text(
  summ_by_period$obs_zero_w_n, summ_by_period$bb_p0_w_n,
  labels = summ_by_period$period, pos = 4, cex = 0.8
)

# 2) Tobit/intreg: zero rate by period
plot(
  summ_by_period$obs_zero_w_n, summ_by_period$tobit_p0_w_n,
  xlab = "Observed zero rate (weighted by n)",
  ylab = "Predicted P(zero) (weighted by n)",
  main = sprintf("Tobit/intreg (cpoint = %s): Zero-rate calibration by period", cpoint),
  pch = 19
)
add_subtitle(good_line_txt)
abline(0, 1, lty = 2)
text(
  summ_by_period$obs_zero_w_n, summ_by_period$tobit_p0_w_n,
  labels = summ_by_period$period, pos = 4, cex = 0.8
)

# 3) Binomial GLM: zero rate by period
plot(
  summ_by_period$obs_zero_w_n, summ_by_period$bin_p0_w_n,
  xlab = "Observed zero rate (weighted by n)",
  ylab = "Predicted P(zero) (weighted by n)",
  main = "Binomial GLM (logit): Zero-rate calibration by period",
  pch = 19
)
add_subtitle(good_line_txt)
abline(0, 1, lty = 2)
text(
  summ_by_period$obs_zero_w_n, summ_by_period$bin_p0_w_n,
  labels = summ_by_period$period, pos = 4, cex = 0.8
)

# 4) Mean share calibration (pick best available alt model)
# Here we plot Beta-binomial mean share; you can swap to binomial/tobit if you want.
plot(
  summ_by_period$obs_share_w_n, summ_by_period$bb_share_w_n,
  xlab = "Observed mean trauma share (weighted by n)",
  ylab = "Predicted mean trauma share (weighted by n)",
  main = "Beta-binomial (glmmTMB): Mean-share calibration by period",
  pch = 19
)
add_subtitle("Good fit: points close to 45° line; errors not concentrated in specific periods.")
abline(0, 1, lty = 2)
text(
  summ_by_period$obs_share_w_n, summ_by_period$bb_share_w_n,
  labels = summ_by_period$period, pos = 4, cex = 0.8
)

par(op)




# -----------------
  
# ============================================================
# APPENDIX: FINAL MODEL + DIAGNOSTICS + CLUSTER-BOOT INFERENCE
# (direction + significance comparable across Tobit vs beta-binomial)
# ============================================================

# --- choose your "final" model here ---
# Recommended from your outputs: standard beta-binomial (m_bb),
# since ZI beta-binomial did not improve AIC/BIC.
final_choice <- "bb"   # "bb" or "bb_zi"

stopifnot(exists("dat2"), exists("f_bin"))
stopifnot(exists("m_bb") || exists("m_bb_zi"))

m_final <- switch(
  final_choice,
  bb    = m_bb,
  bb_zi = m_bb_zi,
  stop("final_choice must be 'bb' or 'bb_zi'")
)

final_label <- switch(
  final_choice,
  bb    = "Beta-binomial (glmmTMB)",
  bb_zi = "Zero-inflated beta-binomial (glmmTMB)"
)

cat("\n====================\nFINAL MODEL\n====================\n")
cat("Final model:", final_label, "\n")
print(summary(m_final))

# ------------------------------------------------------------
# 1) DHARMa diagnostics (titled)
# ------------------------------------------------------------
if (have_pkg("DHARMa")) {
  library(DHARMa)
  
  add_dharma_title <- function(model_label) {
    mtext(paste0("DHARMa residual diagnostics: ", model_label),
          side = 3, line = 1, cex = 1.0, font = 2)
  }
  
  cat("\n====================\nDHARMa diagnostics (final model)\n====================\n")
  res_final <- simulateResiduals(m_final, plot = FALSE)
  plot(res_final)
  add_dharma_title(final_label)
  
  print(testDispersion(res_final))
  print(testZeroInflation(res_final))
} else {
  cat("\nNOTE: Install DHARMa to run residual diagnostics.\n")
}

# ------------------------------------------------------------
# 2) Cluster bootstrap (by co2LF) for robust SE + p-values
# ------------------------------------------------------------
# Why: glmmTMB + sandwich vcov can be fragile; bootstrap is reliable with ~21 clusters.

cluster_boot_glmmTMB <- function(
    model, data, cluster_var,
    B = 400, seed = 1,
    family_obj = betabinomial(link = "logit"),
    control = glmmTMB::glmmTMBControl(
      optCtrl = list(iter.max = 1e4, eval.max = 1e4)
    ),
    verbose = TRUE
) {
  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  set.seed(seed)
  
  cl <- data[[cluster_var]]
  if (!is.factor(cl)) cl <- as.factor(cl)
  clusters <- levels(cl)
  G <- length(clusters)
  
  # point estimates from fitted model
  b_hat <- glmmTMB::fixef(model)$cond
  coef_names <- names(b_hat)
  P <- length(b_hat)
  
  # preserve factor levels (important for matching coefficient columns)
  co2_levels <- if ("co2" %in% names(data) && is.factor(data$co2)) levels(data$co2) else NULL
  
  boot_mat <- matrix(NA_real_, nrow = B, ncol = P)
  colnames(boot_mat) <- coef_names
  n_ok <- 0L
  
  for (b in seq_len(B)) {
    # resample clusters with replacement
    samp <- sample(clusters, size = G, replace = TRUE)
    
    idx <- unlist(lapply(samp, function(g) which(cl == g)), use.names = FALSE)
    d_b <- data[idx, , drop = FALSE]
    
    # enforce factor levels so coef names align
    if (!is.null(co2_levels) && "co2" %in% names(d_b)) {
      d_b$co2 <- factor(d_b$co2, levels = co2_levels)
    }
    
    fit_b <- try(
      glmmTMB::glmmTMB(
        formula(model),
        data = d_b,
        family = family_obj,
        ziformula = formula(model)$ziformula %||% ~ 0,
        control = control
      ),
      silent = TRUE
    )
    
    if (!inherits(fit_b, "try-error")) {
      bb <- glmmTMB::fixef(fit_b)$cond
      # align to original coefficient set
      if (all(coef_names %in% names(bb))) {
        boot_mat[b, ] <- bb[coef_names]
        n_ok <- n_ok + 1L
      }
    }
    
    if (verbose && (b %% 50 == 0)) {
      cat("bootstrap", b, "/", B, " completed; usable fits:", n_ok, "\n")
    }
  }
  
  boot_ok <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
  if (nrow(boot_ok) < max(50, 0.3 * B)) {
    warning("Many bootstrap refits failed or dropped coefficients. Consider reducing model complexity or B.")
  }
  
  se_boot <- apply(boot_ok, 2, sd)
  
  # bootstrap "two-sided sign" p-value: 2*min(P(beta<=0), P(beta>=0))
  p_boot <- sapply(seq_len(P), function(j) {
    v <- boot_ok[, j]
    2 * min(mean(v <= 0), mean(v >= 0))
  })
  p_boot <- pmin(p_boot, 1)
  
  ci_lo <- apply(boot_ok, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_hi <- apply(boot_ok, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  out <- data.frame(
    Estimate = as.numeric(b_hat),
    Boot_SE  = as.numeric(se_boot),
    CI_2.5   = as.numeric(ci_lo),
    CI_97.5  = as.numeric(ci_hi),
    p_boot   = as.numeric(p_boot),
    row.names = coef_names
  )
  
  list(table = out, boot = boot_ok, n_ok = nrow(boot_ok), B = B)
}

# run bootstrap on FINAL model
if (have_pkg("glmmTMB")) {
  library(glmmTMB)
  
  cat("\n====================\nCluster bootstrap inference (final model)\n====================\n")
  boot_res <- cluster_boot_glmmTMB(
    model = m_final,
    data = dat2,
    cluster_var = "co2LF",
    B = 400,      # adjust if you want (e.g., 800 for final write-up)
    seed = 1,
    family_obj = betabinomial(link = "logit"),
    verbose = TRUE
  )
  
  cat("\nUsable bootstrap refits:", boot_res$n_ok, "out of", boot_res$B, "\n")
  final_coef_table <- boot_res$table[order(boot_res$table$p_boot), , drop = FALSE]
  print(final_coef_table)
  
} else {
  cat("\nNOTE: glmmTMB not available; skipping cluster bootstrap.\n")
}

# ------------------------------------------------------------
# Tobit vs Final model: direction + significance (robust)
#   - Tobit: uses cluster-robust p-values already in out_stata
#   - Final model: uses bootstrap p_boot in boot_res$table
# ------------------------------------------------------------

if (exists("out_stata") && exists("boot_res")) {
  
  # --- robustly find the Tobit p-value column ---
  pcol_tobit <- grep("^Pr|p\\b|p_value|pvalue", names(out_stata), value = TRUE)
  if (length(pcol_tobit) == 0) stop("Couldn't find Tobit p-value column in out_stata. Try names(out_stata).")
  pcol_tobit <- pcol_tobit[1]
  
  tobit_tbl <- data.frame(
    term = rownames(out_stata),
    Estimate_tobit = as.numeric(out_stata[,"Estimate"]),
    p_tobit = as.numeric(out_stata[[pcol_tobit]]),
    row.names = NULL
  )
  
  final_tbl <- boot_res$table
  final_tbl <- data.frame(
    term = rownames(final_tbl),
    Estimate_final = as.numeric(final_tbl[,"Estimate"]),
    p_final = as.numeric(final_tbl[,"p_boot"]),
    row.names = NULL
  )
  
  comp <- merge(tobit_tbl, final_tbl, by = "term", all = FALSE)
  
  comp <- within(comp, {
    sign_tobit <- ifelse(Estimate_tobit > 0, "+", "-")
    sign_final <- ifelse(Estimate_final > 0, "+", "-")
    agree_sign <- sign_tobit == sign_final
    sig_tobit_05 <- p_tobit < 0.05
    sig_final_05 <- p_final < 0.05
  })
  
  cat("\n====================\nTobit vs Final model: sign + p-values\n====================\n")
  print(comp[order(comp$p_final), c(
    "term",
    "Estimate_tobit", "p_tobit",
    "Estimate_final", "p_final",
    "agree_sign", "sig_tobit_05", "sig_final_05"
  )])
  
} else {
  cat("\nNOTE: Need both out_stata and boot_res to build comparison table.\n")
}

cat("\nDONE.\n")
