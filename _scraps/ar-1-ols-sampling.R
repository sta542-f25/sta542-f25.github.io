# Function to simulate AR(1) process
simulate_ar_1 <- function(T, b0, b1, s, m0, s0) {
  y <- numeric(T)
  y[1] <- rnorm(1, m0, s0)
  for (t in 2:T) {
    y[t] <- b0 + b1 * y[t - 1] + rnorm(1, 0, s)
  }
  return(y)
}

# --- Parameters ---
set.seed(8675)   
T     <- 5000    
b0    <- 0    
b1    <- 0.5    # 0.5, 0.99, 1, 1.001
s     <- 1      
m0    <- 0      
s0    <- 0      
R     <- 5000   

# --- Storage ---
estimates <- matrix(NA, nrow = R, ncol = 2) # columns: intercept, slope

# --- Simulation loop ---
for (r in 1:R) {
  y <- simulate_ar_1(T, b0, b1, s, m0, s0)
  # Create lagged regressor
  y_lag <- y[-T]
  y_resp <- y[-1]
  # OLS regression
  fit <- lm(y_resp ~ y_lag)
  estimates[r, ] <- coef(fit)
}

# --- Results ---
slope_estimates <- estimates[, 2]

# Histogram of slope estimates
hist(slope_estimates,
     breaks = "Scott",
     main = "Sampling Distribution of OLS Slope Estimates",
     xlab = "Slope estimate",
     col = "lightblue",
     border = "white",
     freq = FALSE)
abline(v = b1, col = "red", lwd = 2)
legend("topright",
       legend = c(
         bquote(beta[0] == .(b0)),
         bquote(beta[1] == .(b1)),
         bquote(sigma == .(s)),
         bquote(T == .(T))
       ),
       bty = "n")
