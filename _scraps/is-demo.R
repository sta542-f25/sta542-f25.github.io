# ==============================================================================
# helper function
# ==============================================================================

beta_is <- function(Mmin, Mmax, true_par, prop_a, prop_b){
  
  estimate <- numeric(Mmax)
  X <- rbeta(Mmax, prop_a, prop_b)
  for(M in Mmin:Mmax){
    x <- X[1:M]
    w <- dbeta(x, true_par, true_par) / dbeta(x, prop_a, prop_b)
    W <- w / sum(w)
    estimate[M] <- sum(log(x) * W)
  }
  
  return(estimate)
}

beta_is_2 <- function(Nrep, M, true_par, prop_a, prop_b){
  
  estimate <- numeric(Nrep)
  for(i in 1:Nrep){
    x <- rbeta(M, prop_a, prop_b)
    w <- dbeta(x, true_par, true_par) / dbeta(x, prop_a, prop_b)
    W <- w / sum(w)
    estimate[i] <- sum(log(x) * W)
  }
  
  return(estimate)
}

beta_is_plot <- function(true_par, prop_a, prop_b, ymax){
  
  curve(dbeta(x, true_par, true_par), from = 0, to = 1, n = 1000, col = "blue", 
        lwd = 3, ylab = "", ylim = c(0, ymax))
  curve(dbeta(x, prop_a, prop_b), from = 0, to = 1, n = 1000, col = "red", lwd = 3, add = TRUE)
  curve(dbeta(x, true_par, true_par) / dbeta(x, prop_a, prop_b), lty = 2,
        from = 0, to = 1, n = 1000, col = "purple", lwd = 3, add = TRUE)
  legend("topright", c("f", "g", "f / g"), lty = c(1, 1, 2), col = c("blue", "red", "purple"), bty = "n", lwd = 3)
  
}

# ==============================================================================
# settings
# ==============================================================================

set.seed(123456)

true_par <- 5
ymax <- 5
Mmin <- 5
Mmax <- 500
Nrep = 1000
M <- 30

# ==============================================================================
# plot alternatives
# ==============================================================================

beta_is_plot(true_par, true_par, true_par, ymax)
beta_is_plot(true_par, 1, 1, ymax)
beta_is_plot(true_par, 20, 20, ymax)
beta_is_plot(true_par, 2, 5, ymax)
beta_is_plot(true_par, 0.3, 0.3, ymax)

# ==============================================================================
# finite-sample behavior
# ==============================================================================

fs_estimate_ideal <- beta_is_2(Nrep, M, true_par, true_par, true_par)
fs_estimate_unif <- beta_is_2(Nrep, M, true_par, 1, 1)
fs_estimate_narrow <- beta_is_2(Nrep, M, true_par, 20, 20)
fs_estimate_left <- beta_is_2(Nrep, M, true_par, 2, 5)
fs_estimate_tails <- beta_is_2(Nrep, M, true_par, 0.3, 0.3)

hist(fs_estimate_ideal, freq = FALSE, breaks = 50, col = rgb(0, 0, 1, alpha = 0.5), main = "Sampling dist of IS", ylab = "estimates")
hist(fs_estimate_unif, freq = FALSE, breaks = 50, add = TRUE, col = rgb(1, 0, 0, alpha = 0.5))
hist(fs_estimate_narrow, freq = FALSE, breaks = 50, add = TRUE, col = rgb(1, 0, 0, alpha = 0.5))
hist(fs_estimate_left, freq = FALSE, breaks = 50, add = TRUE, col = rgb(1, 0, 0, alpha = 0.5))
hist(fs_estimate_tails, freq = FALSE, breaks = 50, add = TRUE, col = rgb(1, 0, 0, alpha = 0.5))

# ==============================================================================
# long run behavior
# ==============================================================================


estimate_ideal <- beta_is(Mmin, Mmax, true_par, true_par, true_par)
estimate_unif <- beta_is(Mmin, Mmax, true_par, 1, 1)
estimate_narrow <- beta_is(Mmin, Mmax, true_par, 20, 20)
estimate_left <- beta_is(Mmin, Mmax, true_par, 2, 5)
estimate_tails <- beta_is(Mmin, Mmax, true_par, 0.3, 0.3)


plot(Mmin:Mmax, estimate_ideal[Mmin:Mmax], type = "l", ylim = c(-1.2, -0.5))
lines(Mmin:Mmax, estimate_unif[Mmin:Mmax])
lines(Mmin:Mmax, estimate_narrow[Mmin:Mmax])
lines(Mmin:Mmax, estimate_left[Mmin:Mmax])
lines(Mmin:Mmax, estimate_tails[Mmin:Mmax])
