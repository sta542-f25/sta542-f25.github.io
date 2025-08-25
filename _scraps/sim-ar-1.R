b0 = 0
b1 = 1.01
m0 = 0
v0 = 1
s0 = 1

T <- 1000

y <- numeric(T)
y[1] = rnorm(1, m0, v0)
for(t in 2:T){
  y[t] = b0 + b1 * y[t - 1] + rnorm(1, 0, sd = s0)
}


plot(1:T, y, type = "l")