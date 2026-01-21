# --- Packages ---
library(readxl)
library(dplyr)

# --- Load data ---
data_url <- "https://zenodo.org/records/8010025/files/S1_file_combined.xlsx?download=1"
tmp_file <- tempfile(fileext = ".xlsx")
download.file(data_url, destfile = tmp_file, mode = "wb", quiet = TRUE)
df_raw <- read_excel(tmp_file)
unlink(tmp_file)

# --- Build ONE analysis dataset (all objects derive from this) ---
dat <- df_raw %>%
  transmute(
    tr_nop = as.numeric(tr_nop),
    ncases = as.numeric(ncases),
    co2LF  = as.factor(co2LF),
    T0n_   = as.numeric(`T0n_`),
    T2be   = as.numeric(T2be),
    T3bm   = as.numeric(T3bm),
    T4bl   = as.numeric(T4bl),
    T5i_   = as.numeric(`T5i_`)
  ) %>%
  filter(
    !is.na(tr_nop), !is.na(ncases), !is.na(co2LF),
    !is.na(T0n_), !is.na(T2be), !is.na(T3bm), !is.na(T4bl), !is.na(T5i_),
    ncases > 0
  )

stopifnot(nrow(dat) == 82)

# --- Model pieces ---
X <- model.matrix(~ T0n_ + T2be + T3bm + T4bl + T5i_, data = dat)  # includes intercept
K <- ncol(X)
N <- nrow(dat)

y <- ifelse(dat$tr_nop > 0, log(dat$tr_nop), NA_real_)
is_cens <- dat$tr_nop == 0
cpoint <- -3

# Stata-normalized aweights: sum(a) = N
a <- dat$ncases
a <- a * N / sum(a)

cluster <- dat$co2LF
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





