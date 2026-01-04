library(depmixS4)

simulate_hmm <- function(n, mu = c(0, 3), sigma = c(1, 1), 
                         trans = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)) {
  state <- numeric(n)
  y <- numeric(n)
  
  # initial state
  state[1] <- sample(1:2, 1, prob = c(0.5, 0.5))
  y[1] <- rnorm(1, mu[state[1]], sigma[state[1]])
  
  for (t in 2:n) {
    state[t] <- sample(1:2, 1, prob = trans[state[t-1], ])
    y[t] <- rnorm(1, mu[state[t]], sigma[state[t]])
  }
  
  data.frame(y = y, state = state)
}

# ==============================================================================
# Sampling distribution (why bimodal?)
# ==============================================================================

Nsim <- 2000
n <- 100

mu1 <- numeric(Nsim)
mu2 <- numeric(Nsim)

for(m in 1:Nsim){
  mydata <- simulate_hmm(n)
  hmm <- depmix(y ~ 1, data = mydata, nstates = 2, family = gaussian())
  fit_hmm <- fit(hmm)
  pars <- getpars(fit_hmm)
  mu1[m] <- pars[7]
  mu2[m] <- pars[9]
}

par(mfrow = c(1, 1))

plot(mu1, mu2, pch = 19, cex = 0.25)



set.seed(42)
n <- 1000
dat <- simulate_hmm(n, trans = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2, byrow = TRUE))

# --- fit HMM -----------------------------------------------------------------
mod <- depmix(response = y ~ 1, family = gaussian(), nstates = 2, data = dat)
fm <- fit(mod, verbose = FALSE)

# --- get smoothed (posterior) probabilities ----------------------------------
post <- posterior(fm, type = "smoothing")
dat$prob_state1 <- post[, 1]
dat$prob_state2 <- post[, 2]

# --- helper: find contiguous segments where true_state == 2 ------------------
flag <- dat$state == 2
r <- rle(flag)

# compute start/end positions for each run
ends_all <- cumsum(r$lengths)
starts_all <- ends_all - r$lengths + 1

# select only runs where flag == TRUE
if (any(r$values)) {
  starts <- starts_all[r$values]
  ends   <- ends_all[r$values]
} else {
  starts <- integer(0)
  ends   <- integer(0)
}

par(mfrow = c(1, 1))

# --- plot --------------------------------------------------------------------
plot(1:nrow(dat), dat$prob_state2, type = "l", lwd = 2, col = "steelblue",
     xlab = "Time", ylab = "",
     main = "Smoothed State Probabilities with True-State Shading")

# add gray shaded bars for true state = 2
usr <- par("usr")  # c(x1, x2, y1, y2)
if (length(starts) > 0) {
  for (i in seq_along(starts)) {
    # use integer boundaries so the bar covers the whole time index
    rect(xleft = starts[i] - 0.5, xright = ends[i] + 0.5,
         ybottom = usr[3], ytop = usr[4],
         col = rgb(1, 0, 0, 0.4), border = NA)
  }
}


# redraw line on top of shaded regions
lines(1:nrow(dat), dat$prob_state2, lwd = 2, col = "blue")

# optional: horizontal reference line at 0.5
abline(h = 0.5, lty = 2, col = "gray40")




#mydata <- simulate_hmm(100)
#hmm <- depmix(y ~ 1, data = mydata, nstates = 2, family = gaussian())
#fit_hmm <- fit(hmm)
#summary(fit_hmm)
#est_states <- posterior(fit_hmm, type = "viterbi")