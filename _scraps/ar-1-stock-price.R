# ==========================================================
# load packages
# ==========================================================

library(quantmod)
library(scoringRules)
library(extraDistr)
library(tidyverse)
library(ggridges)

# ==========================================================
# helper function
# ==========================================================

# online update of AR(1) posterior

nig_update <- function(y_t, x_t, a, b, m, V) {
  # x_t: column vector (p x 1), m: column vector (p x 1), V: (p x p)
  
  # predictive residual
  r <- as.numeric(y_t - t(x_t) %*% m)
  
  # scalar c = 1 + x' V x
  c_val <- as.numeric(1 + t(x_t) %*% V %*% x_t)
  
  # update mean
  m_new <- m + (V %*% x_t / c_val) * r
  
  # update covariance
  Vx <- V %*% x_t
  V_new <- V - (Vx %*% t(Vx)) / c_val
  
  # update shape and scale
  a_new <- a + 0.5
  b_new <- b + 0.5 * (r^2 / c_val)
  
  return(list(a = a_new, b = b_new, m = m_new, V = V_new))
}

# Simulate AR(1)

simulate_ar_1 <- function(T, b0, b1, s, m0, s0) {
  y <- numeric(T)
  y[1] <- rnorm(1, m0, s0)
  for (t in 2:T) {
    y[t] <- b0 + b1 * y[t - 1] + rnorm(1, 0, s)
  }
  return(y)
}

compute_PI_coverage <- function(y, pred_params, alpha = 0.10){
  # y = vector of actual observations
  # pred_params = matrix with columns: df, location, scale
  # alpha = 1 - nominal coverage (default 0.10 for 90% PI)
  
  Tlen <- length(y)
  PI_results <- matrix(NA, nrow = Tlen, ncol = 2)
  colnames(PI_results) <- c(paste0("width_", 100*(1-alpha)), 
                            paste0("covered_", 100*(1-alpha)))
  
  for(t in 1:Tlen){
    df_t <- pred_params[t, "df"]
    mu_t <- pred_params[t, "location"]
    sigma_t <- pred_params[t, "scale"]
    
    # Skip if any parameters are NA
    if(any(is.na(c(df_t, mu_t, sigma_t)))) next
    
    # Compute lower and upper bounds of the PI
    lower <- qlst(alpha/2, df = df_t, mu = mu_t, sigma = sigma_t)
    upper <- qlst(1 - alpha/2, df = df_t, mu = mu_t, sigma = sigma_t)
    
    # Store interval width
    PI_results[t, 1] <- upper - lower
    
    # Store coverage indicator (0/1)
    PI_results[t, 2] <- as.numeric(y[t] >= lower & y[t] <= upper)
  }
  
  return(PI_results)
}

# ==========================================================
# load packages
# ==========================================================

real <- TRUE

if(real){
  getSymbols("AAPL", from = "2000-01-01", to = "2025-09-08", src = "yahoo")
  dates <- index(AAPL)
  
  full_dates <- seq(from = dates[1], to = dates[length(dates)], by = "days")
  
  # Merge with original xts, filling missing dates with NA
  y_full <- merge(AAPL, xts(order.by = full_dates), all = TRUE)
  
  # Plot
  plot(index(y_full), as.numeric(y_full$AAPL.Close), type = "l", col = "blue",
       xlab = "Day", ylab = "Closing price (USD)",
       main = "Apple stock price (2000-2025)")
  grid()
  
  y <- as.numeric(Cl(AAPL))
  
} else {
  set.seed(8675309)
  T <- 5000
  b0 <- 0
  b1 <- 1
  s <- 2
  m0 <- 0
  s0 <- 1
  y <- simulate_ar_1(T, b0, b1, s, m0, s0)
  
  main_title <- paste0("Simulated AR(1) data: ",
                       "β0 = ", b0, ", ",
                       "β1 = ", b1, ", ",
                       "\u03C3 = ", s)
  
  plot(1:T, y, type = "l", main = main_title,
       xlab = "t", ylab = expression(y[t]))
}

Tlen <- length(y)

# ==========================================================
# pre-allocate storage for forecast distributions
# ==========================================================

pred_params <- matrix(NA, nrow = Tlen, ncol = 3,
                      dimnames = list(NULL, c("df", "location", "scale")))

pred_params_iid <- matrix(NA, nrow = Tlen, ncol = 3)
colnames(pred_params_iid) <- c("df", "location", "scale")

# ==========================================================
# prior
# ==========================================================

p <- 2  
m <- matrix(0, nrow = p, ncol = 1)
V <- diag(1, p) 
a <- 3
b <- 1

# ==========================================================
# process data
# ==========================================================

for (t in 2:Tlen) {
  
  y_past <- y[1:(t-1)]
  n <- length(y_past)
  
  x_t <- matrix(c(1, y[t-1]), nrow = p)
  
  df <- 2 * a
  loc <- as.numeric(t(x_t) %*% m)
  scale <- sqrt((b / a) * (1 + t(x_t) %*% V %*% x_t))
  
  pred_params[t, ] <- c(df, loc, scale)
  
  res <- nig_update(y[t], x_t, a, b, m, V)
  a <- res$a; b <- res$b; m <- res$m; V <- res$V
  
  if (n >= 2) {
    ybar <- mean(y_past)
    s <- sd(y_past)
    
    # predictive parameters
    pred_params_iid[t, "df"] <- n - 1
    pred_params_iid[t, "location"] <- ybar
    pred_params_iid[t, "scale"] <- s * sqrt(1 + 1/n)
  }
}

# ==========================================================
# look at mean squared error of predictive mean
# look at mean absolute error of the predictive median
# ==========================================================

mean((y - pred_params[,"location"])^2, na.rm = TRUE)
mean(abs(y - pred_params[,"location"]), na.rm = TRUE)

mean((y - pred_params_iid[,"location"])^2, na.rm = TRUE)
mean(abs(y - pred_params_iid[,"location"]), na.rm = TRUE)

# ==========================================================
# look at calibration
# ==========================================================

PI_matrix <- compute_PI_coverage(y, pred_params, alpha = 0.10)

PI_iid_matrix <- compute_PI_coverage(y, pred_params_iid, alpha = 0.10)

# Average width and empirical coverage
mean(PI_matrix[,1], na.rm = TRUE)         # average interval width
mean(PI_matrix[,2], na.rm = TRUE) 

mean(PI_iid_matrix[,1], na.rm = TRUE)         # average interval width
mean(PI_iid_matrix[,2], na.rm = TRUE) 

# ==========================================================
# look at calibration
# ==========================================================

PIT <- mapply(function(y_t, df_t, mu_t, sigma_t) {
                  if(any(is.na(c(y_t, df_t, mu_t, sigma_t)))) return(NA)
                  plst(q = y_t, df = df_t, mu = mu_t, sigma = sigma_t)
                },
                y, 
                pred_params[,"df"], 
                pred_params[,"location"], 
                pred_params[,"scale"])

hist(PIT, breaks = "Scott", freq = FALSE, col = "lightblue", border = "white")
abline(h = 1)

PIT_iid <- mapply(function(y_t, df_t, mu_t, sigma_t) {
                if(any(is.na(c(y_t, df_t, mu_t, sigma_t)))) return(NA)
                plst(q = y_t, df = df_t, mu = mu_t, sigma = sigma_t)
              },
              y, 
              pred_params_iid[,"df"], 
              pred_params_iid[,"location"], 
              pred_params_iid[,"scale"])

hist(PIT_iid, breaks = "Scott", freq = FALSE, col = "lightblue", border = "white")
abline(h = 1)

# ==========================================================
# look at sharpness
# ==========================================================

hist(pred_params[,"scale"], breaks = "Scott", freq = FALSE,
     main = "Scale parameters of predictive densities",
     xlab = "Student's t scale parameter",
     col = "lightblue", border = "white")

hist(pred_params_iid[,"scale"], breaks = "Scott", freq = FALSE,
     main = "Scale parameters of iid predictive densities",
     xlab = "Student's t scale parameter",
     col = "lightblue", border = "white")


# ==========================================================
# Look at log predictive score
# ==========================================================

# ==========================================================
# Look at continuous ranked probability score
# ==========================================================

# ==========================================================
# waterfall plot of the sequence of densities
# ==========================================================

n_display <- 60
height_scale <- 0.9

valid_rows <- which(!is.na(pred_params[,"df"]) & !is.na(pred_params[,"location"]) & !is.na(pred_params[,"scale"]))

idx <- round(seq(min(valid_rows), max(valid_rows), length.out = n_display))

# --------- build a common x-grid that covers all selected predictive distributions ----------
locs  <- pred_params[idx, "location"]
scales <- pred_params[idx, "scale"]
dfs   <- pred_params[idx, "df"]

x_min <- min(locs - 6 * scales, na.rm = TRUE)
x_max <- max(locs + 6 * scales, na.rm = TRUE)
x_grid <- seq(x_min, x_max, length.out = 10000)   # 300 points is usually smooth enough

# --------- compute densities for each chosen predictive distribution ----------
dens_list <- lapply(seq_along(idx), function(i) {
  ti <- idx[i]
  df  <- as.numeric(pred_params[ti, "df"])
  loc <- as.numeric(pred_params[ti, "location"])
  sc  <- as.numeric(pred_params[ti, "scale"])
  dens <- dlst(x_grid, df, mu = loc, sigma = sc)
  #dens <- dt_locscale(x_grid, df = df, loc = loc, scale = sc)
  data.frame(x = x_grid, ypos = i, density = dens, time_index = ti, loc = loc, scale = sc)
})

dens_df <- do.call(rbind, dens_list) |>
  group_by(time_index) |>
  mutate(density_norm = density / max(density)) |>
  ungroup() |> 
  mutate(height = density_norm * height_scale)

y_positions <- seq_along(idx)
label_vec <- idx  

ggplot(dens_df, aes(x = x, y = ypos, height = height, group = time_index, fill = ypos)) +
  geom_ridgeline(scale = 1, colour = "black", alpha = 0.8, show.legend = FALSE) +
  scale_y_continuous(breaks = y_positions, labels = label_vec, expand = c(0.01, 0)) +
  scale_fill_viridis_c(option = "C") +
  labs(x = "Price (USD)", y = "Period", 
       title = "Waterfall plot of one-step predictive densities",
       subtitle = "Non-standard Student's t from conjugate Bayesian analysis of Gaussian AR(1)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank())

