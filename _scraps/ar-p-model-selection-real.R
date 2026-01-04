library(MASS) # for ginv if needed
library(ggplot2)
library(tidyr)
library(dplyr)

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
    max(Mod(eigen(C)$values)) < 1
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



#---- Posterior predictive density (Student-t) ----
blr_pred_density <- function(y_new, x_new, mu_n, Lambda_n, a_n, b_n) {
  # predictive is Student-t with df = 2 a_n, mean = x' mu_n, scale = sqrt( b_n / a_n * (1 + x' Lambda_n^{-1} x) )
  df <- 2 * a_n
  mean <- as.numeric(t(x_new) %*% mu_n)
  cov_inv <- solve(Lambda_n)
  scale2 <- b_n / a_n * (1 + as.numeric(t(x_new) %*% cov_inv %*% x_new))
  d <- dt((y_new - mean) / sqrt(scale2), df = df, log = TRUE) - 0.5 * log(scale2)
  return(as.numeric(d))
}

#---- Fit AR(p) with Bayesian conjugate prior ----
fit_ar <- function(y, p, mu0, Lambda0, a0, b0) {
  n <- length(y)
  if (n <= p) stop("Not enough data for lag order.")
  
  Y <- y[(p+1):n]
  X <- embed(y, p+1)[, -1, drop = FALSE]  # lagged predictors
  X <- cbind(1, X)  # intercept
  n_eff <- length(Y)
  
  # OLS quantities
  fit <- lm.fit(X, Y)
  sigma2_hat <- sum(fit$residuals^2) / n_eff
  loglik <- -0.5 * n_eff * (log(2*pi*sigma2_hat) + 1)
  k <- ncol(X)
  AIC <- -2*loglik + 2*k
  BIC <- -2*loglik + log(n_eff)*k
  
  # Bayesian update
  post <- blr_nig_update(Y, X, mu0, Lambda0, a0, b0)
  ml <- blr_marginal_likelihood(Y, X, mu0, Lambda0, a0, b0,
                                post$mu_n, post$Lambda_n, post$a_n, post$b_n)
  
  list(AIC = AIC, BIC = BIC, log_marglik = ml$log_marginal_lik,
       post = post, X = X, Y = Y)
}

#---- LFO-CV scoring ----
lfo_cv <- function(y, p, mu0, Lambda0, a0, b0, horizon = 1, type = "point") {
  n <- length(y)
  scores <- c()
  
  for (t0 in (p+1):(n - horizon)) {
    Y_train <- y[1:t0]
    fit <- fit_ar(Y_train, p, mu0, Lambda0, a0, b0)
    post <- fit$post
    
    x_new <- rev(y[(t0 - (p-1)):t0])  # predictor vector
    x_new <- c(1, x_new)
    
    if (type == "point" && horizon == 1) {
      # posterior mean point forecast
      yhat <- as.numeric(t(x_new) %*% post$mu_n)
      err <- (y[t0+1] - yhat)^2
      scores <- c(scores, err)
    }
    if (type == "point" && horizon == 2 && t0+2 <= n) {
      # naive 2-step: plug forecast from step 1 into lag
      yhat1 <- as.numeric(t(x_new) %*% post$mu_n)
      x_new2 <- c(1, yhat1, head(x_new[-1], -1)) # shift lags forward
      yhat2 <- as.numeric(t(x_new2) %*% post$mu_n)
      err <- (y[t0+2] - yhat2)^2
      scores <- c(scores, err)
    }
    if (type == "density" && horizon == 1) {
      lp <- blr_pred_density(y[t0+1], x_new, post$mu_n, post$Lambda_n, post$a_n, post$b_n)
      scores <- c(scores, -lp) # negative log predictive score
    }
  }
  
  mean(scores, na.rm = TRUE)
}

compare_order_selection <- function(T = 200, true_p = 2,
                                    b0 = 0, b = c(0.5, -0.3), s = 1,
                                    m0 = 0, s0 = 1,
                                    max_p = 5) {
  # simulate series
  y <- simulate_ar_p(T, b0, b, s, m0, s0)
  
  # store results
  crit_list <- list()
  logmls <- numeric(max_p)
  
  for (p in 1:max_p) {
    # prior hyperparams: zero mean, identity precision
    mu0 <- matrix(0, p+1, 1)         # intercept + p lags
    Lambda0 <- diag(1, p+1)
    a0 <- 1; b0 <- 1
    
    fit <- fit_ar(y, p, mu0, Lambda0, a0, b0)
    
    # LFO-CV scores
    cv1 <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                  horizon = 1, type = "point")
    cv2 <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                  horizon = 2, type = "point")
    cvdens <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                     horizon = 1, type = "density")
    
    crit_list[[p]] <- c(AIC = fit$AIC,
                        BIC = fit$BIC,
                        log_marglik = fit$log_marglik,
                        CV1 = cv1,
                        CV2 = cv2,
                        CVdens = cvdens)
    logmls[p] <- fit$log_marglik
  }
  
  # posterior model probs (normalize marginal likelihoods)
  logmls_shift <- logmls - max(logmls) # stability
  post_probs <- exp(logmls_shift) / sum(exp(logmls_shift))
  
  # assemble data frame
  df <- as.data.frame(do.call(rbind, crit_list))
  df$Lag <- 1:max_p
  df$PostProb <- post_probs
  df <- df[, c("Lag", "AIC", "BIC", "log_marglik",
               "PostProb", "CV1", "CV2", "CVdens")]
  rownames(df) <- NULL
  
  return(df)
}


compare_order_selection_2 <- function(y, max_p = 5) {
  # simulate series
  T <- length(y)
  
  # store results
  crit_list <- list()
  logmls <- numeric(max_p)
  
  for (p in 1:max_p) {
    # prior hyperparams: zero mean, identity precision
    mu0 <- matrix(0, p+1, 1)         # intercept + p lags
    Lambda0 <- diag(1, p+1)
    a0 <- 1; b0 <- 1
    
    fit <- fit_ar(y, p, mu0, Lambda0, a0, b0)
    
    # LFO-CV scores
    cv1 <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                  horizon = 1, type = "point")
    cv2 <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                  horizon = 2, type = "point")
    cvdens <- lfo_cv(y, p, mu0, Lambda0, a0, b0,
                     horizon = 1, type = "density")
    
    crit_list[[p]] <- c(AIC = fit$AIC,
                        BIC = fit$BIC,
                        log_marglik = fit$log_marglik,
                        CV1 = cv1,
                        CV2 = cv2,
                        CVdens = cvdens)
    logmls[p] <- fit$log_marglik
  }
  
  # posterior model probs (normalize marginal likelihoods)
  logmls_shift <- logmls - max(logmls) # stability
  post_probs <- exp(logmls_shift) / sum(exp(logmls_shift))
  
  # assemble data frame
  df <- as.data.frame(do.call(rbind, crit_list))
  df$Lag <- 1:max_p
  df$PostProb <- post_probs
  df <- df[, c("Lag", "AIC", "BIC", "log_marglik",
               "PostProb", "CV1", "CV2", "CVdens")]
  rownames(df) <- NULL
  
  return(df)
}

# ==============================================================================
# real data
# ==============================================================================

library(quantmod)
getSymbols("AAPL", from = "2023-01-01", to = "2025-09-19", src = "yahoo")
stonks <- as.numeric(Cl(AAPL))
compare_order_selection_2(diff(stonks), 15) |>
  select(Lag, AIC, CV2, BIC, CV1, CVdens, PostProb) |>
  pivot_longer(
    cols = -Lag,
    names_to = "criterion",
    values_to = "value"
  ) |>
  mutate(
    criterion = factor(criterion, 
                       levels = c("AIC", "CV2", "BIC", "CV1", "CVdens", "PostProb"))
  ) |>
  ggplot(aes(x = Lag, y = value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ criterion, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Lag order",
    y = "Criterion value",
    title = "Model selection criteria by lag order"
  )

