#https://guyabel.github.io/fanplot/articles/02_boe.html

library(fanplot)
library(animation)
head(boe)

hard.start = 2005
soft.start = 2010
hard.end = soft.start + 3
soft.end = hard.end - 1

saveGIF(
{

for(y0 in seq(soft.start, soft.end, 0.25)){

  # select relevant data
  #y0 <- 2008
  boe0 <- subset(boe, time0==y0)
  k <- nrow(boe0)
  
  # guess work to set percentiles the BOE are plotting
  p <- seq(0.05, 0.95, 0.05)
  p <- c(0.01, p, 0.99)
  
  # quantiles of split-normal distribution for each probability
  # (row) at each future time point (column)
  cpival <- matrix(NA, nrow = length(p), ncol = k)
  for (i in 1:k)
    cpival[, i] <- qsplitnorm(p, mode = boe0$mode[i], sd = boe0$uncertainty[i], skew = boe0$skew[i])
  
  # past data
  plot(cpi, type = "l", col = "tomato", lwd = 2,
       xlim = c(hard.start, hard.end), ylim = c(-2, 7),
       xaxt = "n", yaxt = "n", ylab="",
       xaxs = "i",
       main = "Bank of England Consumer Price Inflation Forecast")
  
  # background
  rect(y0 - 0.25, par("usr")[3] - 1, hard.end, par("usr")[4] + 1, border = "gray90", col = "gray90")
  
  # add fan
  fan(data = cpival, data.type = "values", probs = p,
      start = y0, frequency = 4,
      anchor = cpi[time(cpi) == y0 - 0.25],
      fan.col = colorRampPalette(c("tomato", "gray90")),
      ln = NULL, rlab = NULL)
  
  # boe aesthetics
  axis(2, at = -2:7, las = 2, tcl = 0.5, labels = FALSE)
  axis(4, at = -2:7, las = 2, tcl = 0.5)
  axis(1, at = hard.start:hard.end, tcl = 0.5)
  axis(1, at = seq(hard.start, hard.end, 0.25), labels = FALSE, tcl = 0.2)
  abline(h = 2)  #boe cpi target
  #abline(v = y0 + 1.75, lty = 2)  #2 year line
}
  
},
movie.name = paste(getwd(), "/slides/images/fan-charts.gif", sep = ""),
interval = 0.3,
ani.width = 1000,
ani.height = 600,
ani.res = 150
)
