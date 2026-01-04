set.seed(123)

# Method-of-moments estimator for MA(1)
mom_ma1 <- function(y) {
  n <- length(y)
  g0 <- mean(y^2)              # sample autocovariance at lag 0
  g1 <- mean(y[-1] * y[-n])    # sample autocovariance at lag 1
  
  # Estimate theta from quadratic relation: g1/g0 = theta / (1 + theta^2)
  r1 <- g1 / g0
  disc <- 1 - 4 * r1^2
  
  if (disc < 0) return(c(theta = NA, sigma2 = NA))  # no real solution
  
  theta1 <- (1 - sqrt(disc)) / (2 * r1)
  theta2 <- (1 + sqrt(disc)) / (2 * r1)
  
  # Choose the invertible root (|theta| < 1)
  theta_hat <- ifelse(abs(theta1) < 1, theta1, theta2)
  
  sigma2_hat <- g0 / (1 + theta_hat^2)
  return(c(theta = theta_hat, sigma2 = sigma2_hat))
}

# Simulation settings
theta_true <- 1.1
sigma_true <- 1
n_reps <- 500
sample_sizes <- c(50, 100, 200, 500, 1000, 2000, 4000, 8000, 16000, 32000)

results <- data.frame(
  n = integer(),
  theta = numeric(),
  sigma2 = numeric()
)

for (n in sample_sizes) {
  for (r in 1:n_reps) {
    # simulate MA(1)
    eps <- rnorm(n + 1, 0, sigma_true)
    y <- eps[2:(n+1)] + theta_true * eps[1:n]
    
    est <- mom_ma1(y)
    results <- rbind(results, data.frame(n = n, theta = est["theta"], sigma2 = est["sigma2"]))
  }
}

# Plotting
par(mfrow = c(1, 2))

boxplot(theta ~ factor(n), data = results,
        main = expression(hat(theta)[1]),
        xlab = "Sample size", ylab = "Estimate",
        col = "lightblue")
abline(h = theta_true, col = "red", lwd = 2)

boxplot(sigma2 ~ factor(n), data = results,
        main = expression(hat(sigma^2)),
        xlab = "Sample size", ylab = "Estimate",
        col = "lightgreen")
abline(h = sigma_true^2, col = "red", lwd = 2)
