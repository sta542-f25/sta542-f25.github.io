simulate_ar_1 <- function(T, b0, b1, s, m0, s0){
  y <- numeric(T)
  y[1] <- rnorm(1, m0, s0)
  for(t in 2:T){
    y[t] <- b0 + b1 * y[t - 1] + rnorm(1, 0, s)
  }
  return(y)
}

ar_1_mean <- function(t, b0, b1, m0){
  if(t == 0){
    return(m0)
  }else{
    return(b0 * sum(b1 ^ (0:(t-1))) + m0 * (b1^t)) 
  }
}

ar_1_var <- function(t, b1, s, s0){
  if(t == 0){
    return(s0^2)
  }else{
    return((s0^2) * (b1^(2*t)) + (s^2) * sum(b1 ^ (2*(0:(t-1)))))
  }
}

ar_1_sd <- function(t, b1, s, s0){
  sqrt(ar_1_var(t, b1, s, s0))
}

b0 = 0
b1 = -1
s = 1
m0 = 0
s0 = 1

range = 0:100
alpha = c(0.01, seq(0.1, 0.9, by = 0.1))

middle <- sapply(range, ar_1_mean, b0, b1, m0)
sds <- sapply(range, ar_1_sd, b1, s, s0)


plot(range, middle, type = "l",
     xaxt = "n", 
     yaxt = "n",
     xlab = "t",
     ylab = expression(y[t]),
     ylim = c(-10, 10), bty = "n",
     col = "white")

for(a in alpha){
  
  U = qnorm(1 - a / 2, mean = middle, sd = sds)
  L = qnorm(a / 2, mean = middle, sd = sds)
  
  polygon(
    c(range, rev(range)),
    c(U, rev(L)),
    col = rgb(1, 0, 0, 0.15),
    border = NA
    )
}

inc = 20
axis(1, pos = 0, at = seq(0, max(range), by = inc), 
     labels = c(NA, seq(inc, max(range), by = inc)))
axis(2, pos = 0)

lines(range, simulate_ar_1(max(range) + 1, b0, b1, s, m0, s0), col = "black", lwd = 2)



#lines(range, middle - sds)
#lines(range, middle + sds)




