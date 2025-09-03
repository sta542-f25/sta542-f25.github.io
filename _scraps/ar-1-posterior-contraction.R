# Function to simulate AR(1)
simulate_ar_1 <- function(T, b0, b1, s, m0, s0) {
  y <- numeric(T)
  y[1] <- rnorm(1, m0, s0)
  for (t in 2:T) {
    y[t] <- b0 + b1 * y[t - 1] + rnorm(1, 0, s)
  }
  return(y)
}

# --- True parameters ---
set.seed(123)
T     <- 200
b0    <- 0
b1    <- 1
s     <- 1
m0    <- 0
s0    <- 0

# --- Simulate one dataset ---
y <- simulate_ar_1(T, b0, b1, s, m0, s0)
y_resp <- y[-1]
X <- cbind(1, y[-T])   # design matrix (intercept + lag)

n <- length(y_resp)

# --- Prior hyperparameters ---
mu0 <- c(0, 0)                  # prior mean
Lambda0 <- diag(2) * 0.01       # prior precision (weakly informative)
a0 <- 2
b0_prior <- 2

# --- Posterior updates ---
XtX <- t(X) %*% X
Xty <- t(X) %*% y_resp

Lambda_n <- Lambda0 + XtX
mu_n <- solve(Lambda_n, Lambda0 %*% mu0 + Xty)
a_n <- a0 + n / 2
b_n <- as.numeric(b0_prior + 0.5 * (t(y_resp) %*% y_resp +
                                      t(mu0) %*% Lambda0 %*% mu0 - t(mu_n) %*% Lambda_n %*% mu_n))

# --- Posterior sampling ---
R <- 5000
beta_draws <- matrix(NA, nrow = R, ncol = 2)

for (r in 1:R) {
  # sample sigma^2 ~ Inv-Gamma(a_n, b_n)
  sigma2 <- 1 / rgamma(1, shape = a_n, rate = b_n)
  # sample beta ~ N(mu_n, sigma^2 * Lambda_n^{-1})
  beta <- mvtnorm::rmvnorm(1, mean = mu_n, sigma = sigma2 * solve(Lambda_n))
  beta_draws[r, ] <- beta
}

# --- Plot posterior of slope ---
hist(beta_draws[, 2],
     breaks = "Scott",
     main = "Posterior of AR(1) slope",
     xlab = "Slope (b1)",
     col = "lightblue", border = "white",
     freq = FALSE,
     xlim = c(-1.5, 1.5))
abline(v = b1, col = "red", lwd = 2)
