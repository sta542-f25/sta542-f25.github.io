# ==============================================================================
# packages
# ==============================================================================

library(dlm)
library(KFAS)

# ==============================================================================
# fake data
# ==============================================================================

T <- 200
true_var <- 1

y <- sin((1:T) / 4) + rnorm(T, sd = sqrt(true_var))


# ==============================================================================
# run the Kalman filter and smoother in the local-level model
# ==============================================================================

# Build system matrices
Zt <- matrix(1)
Tt <- matrix(1)
Rt <- matrix(1)
Qt <- matrix(0.0001)
Ht <- matrix(1)
a1 <- matrix(0)
P1 <- matrix(1)

# Define model
model <- SSModel(y ~ -1 + 
                   SSMcustom(Z = Zt,
                             T = Tt,
                             R = Rt,
                             Q = Qt,
                             a1 = a1,
                             P1 = P1),
                 H = Ht)

# Kalman filter and smoother
ffbs <- KFS(model)

# ==============================================================================
# run the Kalman simulation smoother
# ==============================================================================

state_draws <- simulateSSM(model, type = "states", nsim = 1000)

# ==============================================================================
# plot filtered and smoothed trend
# ==============================================================================

plot(1:T, y, type = "l", lwd = 2, xlab = "Time", ylab = "Value",
     main = "Kalman Filter and Smoother", col = "black")
lines(1:T, c(ffbs$att), col = "red", lwd = 2)
lines(1:T, c(ffbs$alphahat), col = "blue", lwd = 2)


# ==============================================================================
# use dlm package to Gibbs sample in local level
# ==============================================================================


m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1
ndraw = 1000
burn  = 0
thin  = 10

output = dlmGibbsDIG(y,
            dlmModPoly(order = 1, m0 = m0, C0 = C0),
            shape.y = shape.y,
            rate.y = rate.y,
            shape.theta = shape.theta,
            rate.theta = rate.theta,
            n.sample = ndraw,
            thin = thin,
            ind = 1,
            save.states = TRUE)

my_trend <- rowMeans(output$theta[2:(T+1), 1, ])

plot(1:T, y, type = "l", lwd = 2, xlab = "Time", ylab = "Value",
     main = "Kalman Filter and Smoother", col = "black")
#lines(1:T, c(ffbs$att), col = "red", lwd = 2)
#lines(1:T, c(ffbs$alphahat), col = "blue", lwd = 2)
lines(1:T, c(my_trend), col = "red", lwd = 2)
lines(1:T, sin((1:T) / 4), col = "purple", lwd = 3)
