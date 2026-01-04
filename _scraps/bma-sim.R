library(ggplot2)
library(extraDistr)

simulate_ar_p <- function(T, b0, b, s, m0, s0) {
  p <- length(b)              # lag order
  y <- numeric(T + p)         # include space for initial conditions
  
  # simulate initial conditions
  y[1:p] <- rnorm(p, m0, s0)
  
  # simulate forward
  for (t in (p + 1):(T + p)) {
    y[t] <- b0 + sum(b * y[(t - 1):(t - p)]) + rnorm(1, 0, s)
  }
  
  # return only the last T values (exclude initial conditions)
  return(y[(p + 1):(T + p)])
}

simulate_stationary_ar <- function(p, max_attempts = 1000) {
  stopifnot(p >= 1)
  
  is_stationary <- function(b) {
    # companion matrix
    if (length(b) == 1) {
      return(abs(b) < 1)
    }
    C <- rbind(b, cbind(diag(p-1), rep(0, p-1)))
    max(abs(eigen(C)$values)) < 1
  }
  
  for (i in 1:max_attempts) {
    # generate random coefficients, e.g., uniform in [-1,1]
    b_try <- runif(p, -1, 1)
    if (is_stationary(b_try)) return(b_try)
  }
  stop("Failed to generate a stationary AR(p) vector after max_attempts")
}


blr_nig_update <- function(y, X, mu0, Lambda0, a0, b0) {
  y <- as.numeric(y)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  n <- length(y)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  yty <- as.numeric(crossprod(y))
  
  Lambda_n <- XtX + Lambda0
  rhs <- Xty + Lambda0 %*% mu0
  mu_n <- solve(Lambda_n, rhs)
  
  a_n <- a0 + n / 2
  quad0 <- as.numeric(t(mu0) %*% Lambda0 %*% mu0)
  quadn <- as.numeric(t(mu_n) %*% Lambda_n %*% mu_n)
  b_n <- b0 + 0.5 * (yty + quad0 - quadn)
  
  list(mu_n = mu_n,
       Lambda_n = Lambda_n,
       a_n = a_n,
       b_n = b_n)
}


blr_marginal_likelihood <- function(y, X, mu0, Lambda0, a0, b0, 
                                    mu_n, Lambda_n, a_n, b_n) {
  n <- length(y)
  k <- ncol(X)
  
  # determinants via Cholesky for stability
  det_ratio <- determinant(Lambda0, logarithm = TRUE)$modulus -
    determinant(Lambda_n, logarithm = TRUE)$modulus
  
  log_ml <- lgamma(a_n) - lgamma(a0) +
    a0 * log(b0) - a_n * log(b_n) +
    0.5 * det_ratio -
    (n / 2) * log(2 * pi)
  
  ml <- exp(log_ml)
  return(list(log_marginal_lik = as.numeric(log_ml),
              marginal_lik = as.numeric(ml)))
}

set.seed(123)

# --- 1. Simulate AR(5) ---
T <- 50
p_true <- 3
b0_true <- 0
b_true <- simulate_stationary_ar(p_true)
#b_true[1] <- 0.001
s <- 2
m0 <- 0
s0 <- 1

y <- as.numeric(arima.sim(n = T, list(ma = c(0.99, -0.5))))#simulate_ar_p(T, b0_true, b_true, s, m0, s0)

# --- 2. Fit Bayesian regression for p = 1..10 ---
max_p <- 4
mu0_fn <- function(p) matrix(0, p+1, 1)
Lambda0_fn <- function(p) diag(1, p+1)
a0 <- 1
b0_prior <- 1

logmls <- numeric(max_p)
posterior_list <- list()
ols_list <- list()

for (p in 1:max_p) {
  # construct lagged X matrix
  X <- embed(y, p+1)[, -1, drop = FALSE]
  X <- cbind(1, X)
  Y <- y[(p+1):length(y)]
  
  mu0_p <- mu0_fn(p)
  Lambda0_p <- Lambda0_fn(p)
  
  # Bayesian posterior via batch helper
  post <- blr_nig_update(Y, X, mu0_p, Lambda0_p, a0, b0_prior)
  
  # marginal likelihood using your blr_marginal_likelihood
  ml <- blr_marginal_likelihood(Y, X, mu0_p, Lambda0_p, a0, b0_prior,
                                post$mu_n, post$Lambda_n, post$a_n, post$b_n)
  
  logmls[p] <- ml$log_marginal_lik
  posterior_list[[p]] <- post
}

# --- 3. Posterior model probabilities ---
logml_shift <- logmls - max(logmls)
post_probs <- exp(logml_shift) / sum(exp(logml_shift))
best_p <- which.max(post_probs)
best_post <- posterior_list[[best_p]]

# Grid of values to evaluate densities
z_grid <- seq(min(y) - 1, max(y) + 1, length.out = 200)

x_new <- rev(y[(T-max_p):(T-1)])
x_new <- c(1, x_new)
x_new <- matrix(x_new, ncol=1)

x_best <- x_new[1:(best_p+1), , drop=FALSE]

# --- Plug-in density (normal with OLS/best coefficients) ---
mean_plugin <- as.numeric(t(x_best) %*% best_post$mu_n)
sd_plugin <- sqrt(best_post$b_n / best_post$a_n)
dens_plugin_vec <- dnorm(z_grid, mean_plugin, sd_plugin)

# --- Posterior predictive density (Student-t) ---
df_post <- 2 * best_post$a_n
mean_post <- as.numeric(t(x_best) %*% best_post$mu_n)
scale_post <- as.numeric((best_post$b_n / best_post$a_n) * (1 + t(x_best) %*% solve(best_post$Lambda_n) %*% x_best))
dens_postpred_vec <- dlst(z_grid, df_post, mu = mean_post, sigma = sqrt(scale_post), log = FALSE)





# --- Model-averaged posterior predictive ---
dens_modelavg_vec <- numeric(length(z_grid))
for (p in 1:max_p) {
  post <- posterior_list[[p]]
  w <- post_probs[p]
  # subvector of x_new for this model
  x_sub <- x_new[1:(p+1), , drop=FALSE]
  df <- 2 * post$a_n
  mean <- as.numeric(t(x_sub) %*% post$mu_n)
  scale <- as.numeric((post$b_n / post$a_n) * (1 + t(x_sub) %*% solve(post$Lambda_n) %*% x_sub))
  dens_modelavg_vec <- dens_modelavg_vec + w * dlst(z_grid, df, mu = mean, sigma = sqrt(scale), log = FALSE)
}


# --- Base R plot ---
plot(z_grid, dens_plugin_vec, type="l", col="blue", lwd=2,
     ylim=range(c(dens_plugin_vec, dens_postpred_vec, dens_modelavg_vec)),
     xlab="y_next", ylab="Density",
     main="Predictive densities for next observation")
lines(z_grid, dens_postpred_vec, col="red", lwd=2)
lines(z_grid, dens_modelavg_vec, col="darkgreen", lwd=2)
legend("topright", legend=c("Plug-in", "Posterior predictive", "Model-averaged"),
       col=c("blue", "red", "darkgreen"), lwd=2)


plot(1:max_p, post_probs)
