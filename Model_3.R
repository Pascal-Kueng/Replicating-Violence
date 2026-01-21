# ============================================================
# Replicate Stata Table 2, Model 3:
# xi: intreg ltr_nop ltr_nop2 Neolithic Early_Bronze_Age Middle_Bronze_Age Late_Bronze_Age Iron_Age Battle_Sites Royal_Interments Urban i.Region
#     [aw=ncases], robust cl(Region_Period_Cluster_Unit)
#
# NOTE: Estimation code is unchanged vs your Model 1 replica.
# Only covariates + i.Region dummies are added to X (and required vars to dat).
# ============================================================

# --- Packages ---
library(readxl)
library(dplyr)

# --- Load data ---
df_raw <- read_excel("S1_file_combined.xlsx")

# --- Build ONE analysis dataset (all objects derive from this) ---
dat <- df_raw %>%
  transmute(
    tr_nop = as.numeric(tr_nop),
    ncases = as.numeric(ncases),
    
    # cluster var (as in your model 1)
    Region_Period_Cluster_Unit  = as.factor(co2LF),
    
    # Model 3 additions (for xi: i.Region)
    Region    = as.factor(co2),
    
    # period indicators (unchanged)
    Neolithic   = as.numeric(`T0n_`),
    Early_Bronze_Age   = as.numeric(T2be),
    Middle_Bronze_Age   = as.numeric(T3bm),
    Late_Bronze_Age   = as.numeric(T4bl),
    Iron_Age   = as.numeric(`T5i_`),
    
    # Model 3 covariates
    Battle_Sites = as.numeric(battle),
    Royal_Interments  = as.numeric(royal),
    Urban  = as.numeric(urban)
  ) %>%
  filter(
    !is.na(tr_nop), !is.na(ncases), !is.na(Region_Period_Cluster_Unit),
    !is.na(Region),
    !is.na(Neolithic), !is.na(Early_Bronze_Age), !is.na(Middle_Bronze_Age), !is.na(Late_Bronze_Age), !is.na(Iron_Age),
    !is.na(Battle_Sites), !is.na(Royal_Interments), !is.na(Urban),
    ncases > 0
  )

stopifnot(nrow(dat) == 82)

# Match Stata xi: default base category for i.Region (typically "ir" if present)
if ("ir" %in% levels(dat$Region)) dat$Region <- relevel(dat$Region, ref = "ir")

# --- Model pieces ---
# Model 3: add Battle_Sites + Royal_Interments + Urban + i.Region
X <- model.matrix(
  ~ Neolithic + Early_Bronze_Age + Middle_Bronze_Age + Late_Bronze_Age + Iron_Age + Battle_Sites + Royal_Interments + Urban + Region,
  data = dat
)  # includes intercept

K <- ncol(X)
N <- nrow(dat)

y <- ifelse(dat$tr_nop > 0, log(dat$tr_nop), NA_real_)
is_cens <- dat$tr_nop == 0
cpoint <- -3

# Stata-normalized aweights: sum(a) = N
a <- dat$ncases
a <- a * N / sum(a)

cluster <- dat$Region_Period_Cluster_Unit
stopifnot(length(cluster) == N)

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

# --- Per-observation scores U (N x (K+1)) for the *log-likelihood* ---
xb <- as.vector(X %*% beta_hat)
sigma <- sigma_hat

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
# bread = inverse Hessian of NEGATIVE log-likelihood
bread <- solve(fit$hessian)

# meat = sum scores within clusters, crossprod
Sg <- rowsum(U, cluster)   # G x (K+1)
meat <- crossprod(Sg)

V_CR0 <- bread %*% meat %*% bread

# --- Stata-style cluster small-sample correction for ML: multiply by G/(G-1) only ---
G <- nrow(Sg)  # number of clusters
V_stata <- (G/(G - 1)) * V_CR0

# Betas only
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

out_stata

# Optional: report sigma too
c(sigma = sigma_hat, lnsigma = lnsigma_hat)






# =========================
# Post-estimation fit checks
# =========================
# ============================================================
# COMPREHENSIVE MODEL DIAGNOSTICS (for intreg-style left-censor)
# Drop this block AFTER you have:
#   - theta, beta_hat, sigma_hat, xb, y, is_cens, cpoint, a, N, K, G, Vb
# ============================================================

# -------------------------
# 1) Basic fit statistics
# -------------------------
ll_hat <- loglik(theta)

n_unc <- sum(!is_cens)
n_cen <- sum(is_cens)

k_par <- K + 1  # betas + lnsigma
AIC <- -2 * ll_hat + 2 * k_par
BIC <- -2 * ll_hat + log(N) * k_par

# Robust Wald chi2 test of all slopes (exclude intercept)
idx_slope <- which(colnames(X) != "(Intercept)")
b_slope <- beta_hat[idx_slope]
V_slope <- Vb[idx_slope, idx_slope, drop = FALSE]

Wald_chi2 <- as.numeric(t(b_slope) %*% solve(V_slope) %*% b_slope)
df_wald <- length(idx_slope)
p_wald <- pchisq(Wald_chi2, df = df_wald, lower.tail = FALSE)

# Per-observation SD implied by Stata aweights-as-precision in your likelihood
sd_i <- sigma_hat / sqrt(a)

# Predicted probability of being left-censored
p_cens_hat <- pnorm((cpoint - xb) / sd_i)

# Compare predicted vs observed censoring (unweighted and aweight-weighted)
obs_cens_unw  <- mean(is_cens)
pred_cens_unw <- mean(p_cens_hat)

obs_cens_w  <- weighted.mean(is_cens, w = a)
pred_cens_w <- weighted.mean(p_cens_hat, w = a)

# RMSE on uncensored only (descriptive; uses aweights)
res_u <- y[!is_cens] - xb[!is_cens]
rmse_u <- sqrt(weighted.mean(res_u^2, w = a[!is_cens]))

# Pseudo-R2 (McFadden) style: compare to intercept-only intreg on same sample
# (still a likelihood-based summary; not "variance explained")
X0 <- model.matrix(~ 1, data = dat)
K0 <- ncol(X0)

loglik0 <- function(par0) {
  beta0 <- par0[1:K0]
  lns0  <- par0[K0 + 1]
  s0 <- exp(lns0)
  xb0 <- as.vector(X0 %*% beta0)
  
  ll0 <- numeric(N)
  
  idx_u0 <- which(!is_cens)
  r0 <- y[idx_u0] - xb0[idx_u0]
  ll0[idx_u0] <- -0.5 * ( a[idx_u0] * (r0^2) / (s0^2) + log(2*pi) + 2*lns0 - log(a[idx_u0]) )
  
  idx_c0 <- which(is_cens)
  t0 <- (cpoint - xb0[idx_c0]) * sqrt(a[idx_c0]) / s0
  ll0[idx_c0] <- pnorm(t0, log.p = TRUE)
  
  sum(ll0)
}

fit0 <- optim(
  c(0, 0), fn = function(p) -loglik0(p),
  method = "BFGS", hessian = TRUE
)
if (fit0$convergence != 0) stop("Null (intercept-only) optim did not converge")
ll0_hat <- loglik0(fit0$par)

pseudoR2_mcfadden <- 1 - (ll_hat / ll0_hat)

fit_stats <- data.frame(
  N = N,
  uncensored = n_unc,
  left_censored = n_cen,
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
  rmse_uncensored = rmse_u,
  sigma = sigma_hat,
  ll_null = ll0_hat,
  pseudoR2_mcfadden = pseudoR2_mcfadden
)

print(fit_stats)

# -------------------------
# 2) Prediction summaries
# -------------------------

# Expected value of OBSERVED outcome y_obs = max(cpoint, y*)
# where y* ~ N(xb, sd_i^2)
z <- (xb - cpoint) / sd_i
Ey_obs <- pnorm(z) * xb + sd_i * dnorm(z) + (1 - pnorm(z)) * cpoint

pred_summary <- data.frame(
  xb = xb,
  p_cens_hat = p_cens_hat,
  Ey_obs = Ey_obs,
  is_cens = is_cens,
  a = a
)

cat("\nPrediction summary (unweighted):\n")
print(summary(pred_summary[, c("xb", "p_cens_hat", "Ey_obs")]))

cat("\nPrediction summary (aweight-weighted means):\n")
print(c(
  xb = weighted.mean(xb, w = a),
  p_cens_hat = weighted.mean(p_cens_hat, w = a),
  Ey_obs = weighted.mean(Ey_obs, w = a)
))

# -------------------------
# 3) Plots: residuals & fit
# -------------------------

# 3a) Naive uncensored residuals vs fitted (descriptive only in censored models)
plot(xb[!is_cens], res_u,
     xlab = "Fitted xb (latent mean)",
     ylab = "Residual (y - xb) [uncensored only]",
     main = "Uncensored residuals vs fitted")
abline(h = 0, lty = 2)

# 3b) Predicted censoring probability by observed censoring
boxplot(p_cens_hat ~ is_cens,
        names = c("Uncensored", "Censored"),
        ylab = "Predicted P(censored)",
        main = "Predicted censoring probability by observed status")

# 3c) Observed vs predicted censoring probability (jittered)
plot(jitter(as.numeric(is_cens), amount = 0.08), p_cens_hat,
     xaxt = "n",
     xlab = "Observed is_cens (0=uncens, 1=cens)",
     ylab = "Predicted P(censored)",
     main = "Predicted P(censored) vs observed status")
axis(1, at = c(0, 1), labels = c("0", "1"))

# 3d) Calibration of censoring probabilities (deciles)
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

cat("\nCensoring calibration by decile of predicted P(censored):\n")
print(cal)

plot(cal$pred_unw, cal$obs_unw,
     xlab = "Mean predicted P(censored) (bin)",
     ylab = "Observed censoring rate (bin)",
     main = "Censoring calibration (unweighted)")
abline(0, 1, lty = 2)

plot(cal$pred_w, cal$obs_w,
     xlab = "Weighted mean predicted P(censored) (bin)",
     ylab = "Weighted observed censoring rate (bin)",
     main = "Censoring calibration (aweight-weighted)")
abline(0, 1, lty = 2)

# -------------------------
# 4) Distributional diagnostics:
#    Randomized quantile residuals for censored normal models
# -------------------------
# If model is correct, rq ~ approximately N(0,1)
set.seed(1)

t_c <- (cpoint - xb) / sd_i
Phi_c <- pnorm(t_c)

u <- numeric(N)
u[!is_cens] <- pnorm((y[!is_cens] - xb[!is_cens]) / sd_i[!is_cens])
u[is_cens]  <- runif(sum(is_cens), min = 0, max = Phi_c[is_cens])

# Guard against numerical edge cases (exact 0/1)
eps <- 1e-12
u <- pmin(pmax(u, eps), 1 - eps)

rq <- qnorm(u)

qqnorm(rq, main = "Randomized quantile residuals (should be ~N(0,1))")
qqline(rq)

hist(rq, main = "Randomized quantile residuals", xlab = "rq")

plot(xb, rq,
     xlab = "Fitted xb (latent mean)",
     ylab = "Randomized quantile residual (rq)",
     main = "RQ residuals vs fitted")
abline(h = 0, lty = 2)

# -------------------------
# 5) Influence / outliers (simple, descriptive)
# -------------------------

# Standardized residuals for uncensored obs (uses sd_i)
std_res_u <- res_u / sd_i[!is_cens]
plot(xb[!is_cens], std_res_u,
     xlab = "Fitted xb (latent mean)",
     ylab = "Std residual (uncensored)",
     main = "Standardized residuals (uncensored)")
abline(h = c(-2, 0, 2), lty = c(2, 2, 2))

# Identify a few largest absolute standardized residuals among uncensored
ord <- order(abs(std_res_u), decreasing = TRUE)
top_k <- min(5, length(ord))
top_outliers <- data.frame(
  row = which(!is_cens)[ord[1:top_k]],
  xb = xb[!is_cens][ord[1:top_k]],
  y = y[!is_cens][ord[1:top_k]],
  resid = res_u[ord[1:top_k]],
  std_resid = std_res_u[ord[1:top_k]],
  weight_a = a[!is_cens][ord[1:top_k]],
  Region_Period_Cluster_Unit = cluster[which(!is_cens)[ord[1:top_k]]]
)

cat("\nTop uncensored outliers by |standardized residual|:\n")
print(top_outliers)

# -------------------------
# 6) (Optional) Compare to a no-covariate censoring-only baseline
# -------------------------
# A very rough benchmark: how well does xb separate censored vs uncensored?
# Not a classification model, but can show whether censoring is learnable from X.
o <- order(xb)
plot(xb[o], as.numeric(is_cens)[o],
     xlab = "xb (sorted)",
     ylab = "Observed is_cens (0/1)",
     main = "Observed censoring vs xb (sorted)")
lines(xb[o], p_cens_hat[o], lty = 1)  # overlay predicted prob (scale matches y in [0,1])
abline(h = obs_cens_unw, lty = 2)




tapply(a, is_cens, summary)


o <- order(a, decreasing = TRUE)
head(data.frame(a=a[o], xb=xb[o], is_cens=is_cens[o], p_cens=p_cens_hat[o]), 15)




# Which observations contribute most to "surprise" given censoring status?
# For uncensored obs, look at predicted P(uncensored) = 1 - p_cens_hat
p_unc_hat <- 1 - p_cens_hat
surprise <- ifelse(is_cens, -log(p_cens_hat), -log(p_unc_hat))  # NLL contribution of censoring indicator only

# weight it (since your likelihood is weight-sensitive)
w_surprise <- a * surprise

o2 <- order(w_surprise, decreasing = TRUE)
head(data.frame(
  row = o2,
  a = a[o2],
  xb = xb[o2],
  is_cens = is_cens[o2],
  p_cens = p_cens_hat[o2],
  surprise = surprise[o2],
  w_surprise = w_surprise[o2]
), 10)


report::report_packages()

