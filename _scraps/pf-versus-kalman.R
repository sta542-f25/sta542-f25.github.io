library(MASS)
library(mvtnorm)

plot_bvn <- function(mu, Sigma, levels = c(.25, .5, .75, .9, .99),
                     xlim = NULL, ylim = NULL, n = 100, main = "") {
  
  stopifnot(length(mu) == 2,
            all(dim(Sigma) == c(2, 2)))
  
  # Convert probability levels → density levels
  detS <- det(Sigma)
  fmax <- 1 / (2 * pi * sqrt(detS))
  dens_levels <- fmax * exp(-0.5 * qchisq(levels, df = 2))
  
  # Default limits
  if (is.null(xlim)) xlim <- mu[1] + c(-4, 4) * sqrt(Sigma[1,1])
  if (is.null(ylim)) ylim <- mu[2] + c(-4, 4) * sqrt(Sigma[2,2])
  
  # Grid
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  
  # Density function
  invS <- solve(Sigma)
  const <- 1 / (2 * pi * sqrt(detS))
  dens <- function(pt) {
    dif <- pt - mu
    const * exp(-0.5 * t(dif) %*% invS %*% dif)
  }
  
  # Correct z grid using outer
  z <- outer(x, y, Vectorize(function(a, b) dens(c(a, b))))
  
  # Plot
  contour(x, y, z, levels = dens_levels,
          xlab = "x", ylab = "y", main = main)
}



simulate_local_level_bivariate <- function(n, Q, R, x0 = c(0, 0)) {
  # Input checks
  if (!is.matrix(Q) || !all(dim(Q) == c(2,2)))
    stop("Q must be a 2x2 covariance matrix.")
  if (!is.matrix(R) || !all(dim(R) == c(2,2)))
    stop("R must be a 2x2 covariance matrix.")
  if (det(Q) <= 0 || det(R) <= 0)
    stop("Q and R must be positive definite.")
  
  # Precompute Cholesky factors
  C_Q <- chol(Q)
  C_R <- chol(R)
  
  # Allocate
  alpha <- matrix(0, nrow = 2, ncol = n)
  y      <- matrix(0, nrow = 2, ncol = n)
  
  # Initial state
  alpha[, 1] <- x0 + C_Q %*% rnorm(2)
  y[, 1]     <- alpha[, 1] + C_R %*% rnorm(2)
  
  # Recursion
  for (t in 2:n) {
    alpha[, t] <- alpha[, t-1] + C_Q %*% rnorm(2)
    y[, t]     <- alpha[, t]   + C_R %*% rnorm(2)
  }
  
  list(states = alpha, observations = y)
}


pf_step <- function(
    particles,      # matrix: 2 x N, each column a particle state
    weights,        # numeric vector length N, normalized
    z,              # new bearing observation (radians)
    V = diag(2),           # std dev of acceleration noise (m/s^2)
    W = diag(2),              # obs noise SD (radians)
    ess_threshold = 0    #
) {
  
  N <- ncol(particles)

  particles_pred <- matrix(0, 2, N)
  log_w <- numeric(N)
  for(i in 1:N){
    particles_pred[, i] <- particles[, i] + mvrnorm(1, numeric(2), W)
    log_w[i] <- log(weights[i]) + dmvnorm(z, particles_pred[, i], V, log = TRUE)
  }
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

set.seed(12345)

T <- 10

fake_data <- simulate_local_level_bivariate(T, diag(2), diag(2))
y <- t(fake_data$observations)

# Build model matrices
Zt <- diag(2)
Tt <- diag(2)
Rt <- diag(2)
Qt <- diag(2)
Ht <- diag(2)
a1 <- fake_data$states[, 1]
P1 <- diag(2)

# Define model
model <- SSModel(y ~ -1 + 
                   SSMcustom(Z = Zt,
                             T = Tt,
                             R = Rt,
                             Q = Qt,
                             a1 = a1,
                             P1 = P1),
                 H = Ht)

ffbs <- KFS(model)

Npart <- 200
max_cex <- 2

all_particles <- array(0, dim = c(length(a1), Npart, T))
all_weights <- matrix(0, T, Npart)

particles <- matrix(a1, length(a1), Npart)
weights <- rep(1 / Npart, Npart)

all_particles[, , 1] <- particles
all_weights[1, ] <- weights

for(t in 1:T){
  new_stuff <- pf_step(particles, weights, y[t, ], ess_threshold = Npart)
  particles <- new_stuff$particles
  weights <- new_stuff$weights
  esti <- particles %*% weights
  all_particles[, , t] <- particles
  all_weights[t, ] <- weights
  #ranks <- rank(weights)
  plot_bvn(ffbs$att[t, ], ffbs$Ptt[, , t], xlim = c(-5, 5), ylim = c(-5, 5), main = paste("t =", t))
  points(ffbs$att[t, 1], ffbs$att[t, 2], pch = 19, cex = 2, col = "blue")
  points(t(particles), cex = max_cex, bty = "n", col = rgb(1, 0, 0, alpha = 0.25))
  points(t(particles), pch = 19, cex = max_cex, bty = "n", col = rgb(1, 0, 0, alpha = weights))
  points(esti[1], esti[2], col = "red", pch = 4)
  #points(fake_data$states[1, t], fake_data$states[3, t], col = "blue", cex = 2, pch = 4)
}