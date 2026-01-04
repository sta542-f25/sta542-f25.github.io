# Helper: wrap angle to (-pi, pi]
wrapToPi <- function(a) {
  a <- (a + pi) %% (2*pi) - pi
  return(a)
}

# Simulate bearings-only tracking: constant-velocity model with white acceleration noise
simulate_bearings <- function(N,
                              q = 0.001^2,           # std dev of acceleration noise (m/s^2)
                              r = 0.005^2,           # std dev of bearing noise (radians) (default 1 degree)
                              x0 = c(10, 5, 10, 2)  # initial state: [x, vx, y, vy]
                              ) { 
  
  if(length(x0) != 4) stop("x0 must be length 4: c(x, vx, y, vy)")
  
  # State transition matrices for constant velocity in 2D
  F <- matrix(0, 4, 4)
  F[1,1] <- 1; F[1,2] <- 1
  F[2,2] <- 1
  F[3,3] <- 1; F[3,4] <- 1
  F[4,4] <- 1
  
  G <- matrix(0, 4, 2)
  G[1, 1] <- 0.5
  G[2, 1] <- 1
  G[3, 2] <- 0.5
  G[4, 2] <- 1
  
  # Preallocate
  x <- matrix(0, nrow = 4, ncol = N)   # true states
  z <- numeric(N)                      # noisy bearings (radians)
  z_clean <- numeric(N)                # noiseless bearings (for reference)
  
  # initial state
  x[,1] <- x0
  
  # simulate
  for(t in 2:N) {
    x[,t] <- F %*% x[,t-1] + G %*% rnorm(2, mean = 0, sd = sqrt(q))
  }
  # produce measurements for each time
  for(t in 1:N) {
    bearing_clean <- atan2(x[3,t], x[1,t])
    z_clean[t] <- wrapToPi(bearing_clean)
    # additive Gaussian noise on angle (wrapped)
    z[t] <- wrapToPi(bearing_clean + rnorm(1, 0, sqrt(r)))
  }
  
  return(list(
    states = x,              # 4 x N matrix: [x, vx, y, vy]
    bearing = z,             # noisy bearings (radians)
    bearing_clean = z_clean  # noise-free bearings (radians)
  ))
}


set.seed(1)
sim <- simulate_bearings(N = 60)

# plot true trajectory and bearing lines from observer (origin)
x <- sim$states
plot(x[1,], x[3,], type = "b", pch = 20, xlab = "x", ylab = "y", main = "True trajectory and bearings",
     xlim = c(0, 100), ylim = c(0, 100))
points(0,0, pch = 4, cex = 1.5) # observer
# draw bearing lines (noisy measurements) from observer
for(t in 1:ncol(x)) {
  ang <- sim$bearing[t]
  lines(c(0, 1200*cos(ang)), c(0, 1200*sin(ang)), lty = 2, col = "gray")
}
points(x[1,], x[3,], pch = 20, col = "blue")


