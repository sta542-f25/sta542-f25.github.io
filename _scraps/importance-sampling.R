# -----------------------------------------------
# Importance Sampling Demo (Normal vs t)
# -----------------------------------------------

set.seed(123)

# Helper: importance sampling function
importance_sampling <- function(target_dens, proposal_dens, proposal_rvs, h, n = 1e5, normalize = FALSE) {
  x <- proposal_rvs(n)
  w <- target_dens(x) / proposal_dens(x)
  if (normalize) {
    w <- w / sum(w)
    estimate <- sum(w * h(x))
  } else {
    estimate <- mean(w * h(x))
  }
  list(estimate = estimate, ess = (sum(w)^2) / sum(w^2))
}

# Define densities and samplers
dt3 <- function(x) dt(x, df = 3)
rt3 <- function(n) rt(n, df = 3)
dnorm0 <- function(x) dnorm(x)
rnorm0 <- function(n) rnorm(n)

# Target function: E[h(X)] under target
h <- function(x) x^2

# True expectations
true_norm <- integrate(function(x) x^2 * dnorm0(x), -Inf, Inf)$value
true_t3   <- integrate(function(x) x^2 * dt3(x), -Inf, Inf)$value

# -----------------------------------------------
# Case 1: Target = Normal, Proposal = t(3)
# -----------------------------------------------
cat("Case 1: Target = Normal, Proposal = t(3)\n")
res1_std  <- importance_sampling(dnorm0, dt3, rt3, h, normalize = FALSE)
res1_norm <- importance_sampling(dnorm0, dt3, rt3, h, normalize = TRUE)
cat(sprintf("Standard IS: estimate = %.3f, ESS = %.0f\n", res1_std$estimate, res1_std$ess))
cat(sprintf("Self-norm  IS: estimate = %.3f, ESS = %.0f\n\n", res1_norm$estimate, res1_norm$ess))

# -----------------------------------------------
# Case 2: Target = t(3), Proposal = Normal
# -----------------------------------------------
cat("Case 2: Target = t(3), Proposal = Normal\n")
res2_std  <- importance_sampling(dt3, dnorm0, rnorm0, h, normalize = FALSE)
res2_norm <- importance_sampling(dt3, dnorm0, rnorm0, h, normalize = TRUE)
cat(sprintf("Standard IS: estimate = %.3f, ESS = %.0f\n", res2_std$estimate, res2_std$ess))
cat(sprintf("Self-norm  IS: estimate = %.3f, ESS = %.0f\n\n", res2_norm$estimate, res2_norm$ess))

# -----------------------------------------------
# Visualization: Weight distributions
# -----------------------------------------------
par(mfrow = c(1,2))
n <- 1e4
x1 <- rt3(n)
w1 <- dnorm0(x1) / dt3(x1)
hist(log(w1), breaks = 40, main = "Weights: Normal / t(3)", xlab = "log weight")

x2 <- rnorm0(n)
w2 <- dt3(x2) / dnorm0(x2)
hist(log(w2), breaks = 40, main = "Weights: t(3) / Normal", xlab = "log weight")
