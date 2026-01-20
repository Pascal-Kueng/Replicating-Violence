library(readxl)
library(dplyr)

viol_data <- read_excel("/Users/jasmin/Library/Mobile Documents/com~apple~CloudDocs/PhD/Replication challenge/S1_file_combined.xlsx")

# Explanation: Some of the rows in the data set are not merged despite all variables except for ratio/percentage one's being equal
# To test for robustness, we ran model 3 on the original and merged data set and compared results

#####################################
##### Merging of duplicate rows #####
#####################################

# Check for duplicates (ignoring ncases, trauma, tr_nop, rcc)
viol_data[duplicated(viol_data %>% select(-ncases, -trauma, -tr_nop, -rcc)) | 
            duplicated(viol_data %>% select(-ncases, -trauma, -tr_nop, -rcc), fromLast = TRUE), ]
# For two sites (Tell Brak and Titriş Höyük) two rows exist that share all characteristics except for ncases, trauma, tr_nop, and rcc 

# Merge duplicates with weighted averages
viol_data_merged <- viol_data %>%
  group_by(across(-c(ncases, trauma, tr_nop, rcc))) %>%
  summarise(
    trauma = sum(ncases * trauma) / sum(ncases),
    tr_nop = sum(ncases * tr_nop) / sum(ncases),
    rcc = sum(ncases * rcc) / sum(ncases),
    ncases = sum(ncases),
    .groups = 'drop'
  )

###################################
##### Model 3 - original data #####
###################################

X1 <- model.matrix(
  ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2,
  data = viol_data
)
K1 <- ncol(X1)
N1 <- nrow(viol_data)
y1 <- ifelse(viol_data$tr_nop > 0, log(viol_data$tr_nop), NA_real_)
is_cens1 <- viol_data$tr_nop == 0
cpoint1 <- -3

a1 <- viol_data$ncases
a1 <- a1 * N1 / sum(a1)
cluster1 <- viol_data$co2LF
stopifnot(length(cluster1) == N1)

loglik1 <- function(par) {
  beta <- par[1:K1]
  lns  <- par[K1 + 1]
  s <- exp(lns)
  xb <- as.vector(X1 %*% beta)
  
  ll <- numeric(N1)
  
  idx_u <- which(!is_cens1)
  r <- y1[idx_u] - xb[idx_u]
  ll[idx_u] <- -0.5 * ( a1[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a1[idx_u]) )
  
  idx_c <- which(is_cens1)
  t <- (cpoint1 - xb[idx_c]) * sqrt(a1[idx_c]) / s
  ll[idx_c] <- pnorm(t, log.p = TRUE)
  
  sum(ll)
}

start1 <- c(rep(0, K1), 0)
fit1 <- optim(start1, fn = function(p) -loglik1(p), method = "BFGS", hessian = TRUE)
if (fit1$convergence != 0) stop("Model 3 original: optim did not converge")

theta1 <- fit1$par
beta_hat1 <- theta1[1:K1]
sigma_hat1 <- exp(theta1[K1 + 1])

xb1 <- as.vector(X1 %*% beta_hat1)
U1 <- matrix(0, nrow = N1, ncol = K1 + 1)
colnames(U1) <- c(colnames(X1), "lnsigma")

idx_u1 <- which(!is_cens1)
r1 <- y1[idx_u1] - xb1[idx_u1]
U1[idx_u1, 1:K1]   <- (a1[idx_u1] * r1 / sigma_hat1^2) * X1[idx_u1, , drop = FALSE]
U1[idx_u1, K1 + 1] <- -1 + (a1[idx_u1] * r1^2 / sigma_hat1^2)

idx_c1 <- which(is_cens1)
t1 <- (cpoint1 - xb1[idx_c1]) * sqrt(a1[idx_c1]) / sigma_hat1
log_phi1 <- dnorm(t1, log = TRUE)
log_Phi1 <- pnorm(t1, log.p = TRUE)
lambda1 <- exp(pmin(log_phi1 - log_Phi1, 700))
U1[idx_c1, 1:K1]   <- -(lambda1 * sqrt(a1[idx_c1]) / sigma_hat1) * X1[idx_c1, , drop = FALSE]
U1[idx_c1, K1 + 1] <- -(lambda1 * t1)

bread1 <- solve(fit1$hessian)
Sg1 <- rowsum(U1, cluster1)
meat1 <- crossprod(Sg1)
V_CR0_1 <- bread1 %*% meat1 %*% bread1

G1 <- nrow(Sg1)
V_stata1 <- (G1/(G1 - 1)) * V_CR0_1

Vb1 <- V_stata1[1:K1, 1:K1, drop = FALSE]
se1 <- sqrt(diag(Vb1))
z1  <- beta_hat1 / se1
p1  <- 2 * pnorm(-abs(z1))

model3_original <- data.frame(
  Estimate     = beta_hat1,
  `Std. Error` = se1,
  `z value`    = z1,
  `Pr(>|z|)`   = p1,
  row.names    = colnames(X1)
)

#################################
##### Model 3 - merged data #####
#################################

X2 <- model.matrix(
  ~ T0n_ + T2be + T3bm + T4bl + T5i_ + battle + royal + urban + co2,
  data = viol_data_merged
)
K2 <- ncol(X2)
N2 <- nrow(viol_data_merged)
y2 <- ifelse(viol_data_merged$tr_nop > 0, log(viol_data_merged$tr_nop), NA_real_)
is_cens2 <- viol_data_merged$tr_nop == 0
cpoint2 <- -3

a2 <- viol_data_merged$ncases
a2 <- a2 * N2 / sum(a2)
cluster2 <- viol_data_merged$co2LF
stopifnot(length(cluster2) == N2)

loglik2 <- function(par) {
  beta <- par[1:K2]
  lns  <- par[K2 + 1]
  s <- exp(lns)
  xb <- as.vector(X2 %*% beta)
  
  ll <- numeric(N2)
  
  idx_u <- which(!is_cens2)
  r <- y2[idx_u] - xb[idx_u]
  ll[idx_u] <- -0.5 * ( a2[idx_u] * (r^2) / (s^2) + log(2*pi) + 2*lns - log(a2[idx_u]) )
  
  idx_c <- which(is_cens2)
  t <- (cpoint2 - xb[idx_c]) * sqrt(a2[idx_c]) / s
  ll[idx_c] <- pnorm(t, log.p = TRUE)
  
  sum(ll)
}

start2 <- c(rep(0, K2), 0)
fit2 <- optim(start2, fn = function(p) -loglik2(p), method = "BFGS", hessian = TRUE)
if (fit2$convergence != 0) stop("Model 3 merged: optim did not converge")

theta2 <- fit2$par
beta_hat2 <- theta2[1:K2]
sigma_hat2 <- exp(theta2[K2 + 1])

xb2 <- as.vector(X2 %*% beta_hat2)
U2 <- matrix(0, nrow = N2, ncol = K2 + 1)
colnames(U2) <- c(colnames(X2), "lnsigma")

idx_u2 <- which(!is_cens2)
r2 <- y2[idx_u2] - xb2[idx_u2]
U2[idx_u2, 1:K2]   <- (a2[idx_u2] * r2 / sigma_hat2^2) * X2[idx_u2, , drop = FALSE]
U2[idx_u2, K2 + 1] <- -1 + (a2[idx_u2] * r2^2 / sigma_hat2^2)

idx_c2 <- which(is_cens2)
t2 <- (cpoint2 - xb2[idx_c2]) * sqrt(a2[idx_c2]) / sigma_hat2
log_phi2 <- dnorm(t2, log = TRUE)
log_Phi2 <- pnorm(t2, log.p = TRUE)
lambda2 <- exp(pmin(log_phi2 - log_Phi2, 700))
U2[idx_c2, 1:K2]   <- -(lambda2 * sqrt(a2[idx_c2]) / sigma_hat2) * X2[idx_c2, , drop = FALSE]
U2[idx_c2, K2 + 1] <- -(lambda2 * t2)

bread2 <- solve(fit2$hessian)
Sg2 <- rowsum(U2, cluster2)
meat2 <- crossprod(Sg2)
V_CR0_2 <- bread2 %*% meat2 %*% bread2

G2 <- nrow(Sg2)
V_stata2 <- (G2/(G2 - 1)) * V_CR0_2

Vb2 <- V_stata2[1:K2, 1:K2, drop = FALSE]
se2 <- sqrt(diag(Vb2))
z2  <- beta_hat2 / se2
p2  <- 2 * pnorm(-abs(z2))

model3_merged <- data.frame(
  Estimate     = beta_hat2,
  `Std. Error` = se2,
  `z value`    = z2,
  `Pr(>|z|)`   = p2,
  row.names    = colnames(X2)
)

######################################################
##### Comparison of estimates (robustness check) #####
######################################################
model3_original
model3_merged

round(model3_original, 3)
round(model3_merged, 3)

round(model3_original-model3_merged, 2)

# Results do not change qualitatively and estimates only slightly differ
# z-values for early and middle bronze age show the largest differences; all duplicated rows were in early bronze age (2be)