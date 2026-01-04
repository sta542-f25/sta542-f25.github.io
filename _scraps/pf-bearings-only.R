# Helper: wrap angle to (-pi, pi]
wrapToPi <- function(a) {
  a <- (a + pi) %% (2*pi) - pi
  return(a)
}

# Simulate bearings-only tracking: constant-velocity model with white acceleration noise
simulate_bearings <- function(N,
                              q = 0.001^2,           # std dev of acceleration noise (m/s^2)
                              r = 0.005^2,           # std dev of bearing noise (radians) (default 1 degree)
                              x0 = c(-0.05, 0.001, 0.7, -0.055)  # initial state: [x, vx, y, vy]
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

pf_step <- function(
    particles,      # matrix: 4 x N, each column a particle state
    weights,        # numeric vector length N, normalized
    z,              # new bearing observation (radians)
    q = 0.001^2,           # std dev of acceleration noise (m/s^2)
    r = 0.005^2,              # obs noise SD (radians)
    ess_threshold = 0    #
  ) {
  
  N <- ncol(particles)
  
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
  
  # simulate importance distribution
  particles_pred <- F %*% particles + G %*% matrix(rnorm(2*N, mean = 0, sd = sqrt(q)), 2, N)
  
  # compute new weights
  bearing_pred <- atan2(particles_pred[3,], particles_pred[1,])
  
  log_w <- log(weights) + dnorm(rep(z, N), mean = bearing_pred, sd = sqrt(r), log = TRUE)
  log_w <- log_w - max(log_w)   # stabilize
  w <- exp(log_w)
  w <- w / sum(w)
  
  # compute ESS
  
  ess <- 1 / sum(w^2)
  
  # resample?
  
  if(ess < ess_threshold) {
    idx <- sample.int(N, size = N, replace = TRUE, prob = w)
    particles_new <- particles_pred[, idx]
    w_new <- rep(1 / N, N)
    resampled <- TRUE
  } else {
    particles_new <- particles_pred
    w_new <- w
    resampled <- FALSE
  }
  
  ## --- 5. Return --------------------------------------------------------------
  list(
    particles = particles_new,
    weights   = w_new,
    ess       = ess,
    resampled = resampled
  )
}

T <- 10
max_cex <- 1
Npart <- 200
q = 0.1
x0 <- c(-0.05, 0.001, 0.7, -0.055)

all_particles <- array(0, dim = c(length(x0), Npart, T))
all_weights <- matrix(0, T, Npart)

particles <- matrix(x0, length(x0), Npart)
particles[1, ] <- particles[1, ] + 0.05 * rnorm(Npart)
particles[3, ] <- particles[3, ] + 0.05 * rnorm(Npart)
weights <- rep(1 / Npart, Npart)

all_particles[, , 1] <- particles
all_weights[1, ] <- weights

set.seed(123456)

fake_data <- simulate_bearings(T, q = 0.001)
z <- fake_data$bearing

for(t in 2:T){
  new_stuff <- pf_step(particles, weights, z[t], ess_threshold = Npart, q = 0.001)
  particles <- new_stuff$particles
  weights <- new_stuff$weights
  all_particles[, , t] <- particles
  all_weights[t, ] <- weights
  plot(0, 0, xlim = c(-1, 1), ylim = c(-1, 1), col = "red", pch = 19, cex = 2,
       xaxt = "n", yaxt = "n", bty = "n")
  axis(1, pos = 0)
  axis(2, pos = 0)
  points(particles[1, ], particles[3, ], cex = max_cex, bty = "n", col = rgb(1, 0, 0, alpha = 0.25))
  points(particles[1, ], particles[3, ], pch = 19, cex = max_cex, bty = "n", col = rgb(1, 0, 0, alpha = weights))
  points(fake_data$states[1, t], fake_data$states[3, t], col = "blue", cex = 2, pch = 4)
}


