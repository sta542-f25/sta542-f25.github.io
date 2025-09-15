library(shiny)

# ---------- AR(2) Simulation ----------
simulate_ar_2 <- function(T, b0, b1, b2, s, m0, s0){
  y <- numeric(T)
  # independent normal initial conditions
  y[1] <- rnorm(1, m0, s0)
  y[2] <- rnorm(1, m0, s0)
  for(t in 3:T){
    y[t] <- b0 + b1 * y[t-1] + b2 * y[t-2] + rnorm(1, 0, s)
  }
  return(y)
}

# ---------- AR(2) State Space Recursions ----------
ar2_state_space <- function(T, b0, b1, b2, sigma, m0, s0) {
  # Companion matrix
  A <- matrix(c(b1, b2,
                1,  0), nrow = 2, byrow = TRUE)
  
  # Constant vector
  c <- c(b0, 0)
  
  # Noise covariance
  Sigma_u <- matrix(c(sigma^2, 0,
                      0,      0), nrow = 2, byrow = TRUE)
  
  # Initialize mean vector and covariance matrix
  mu <- matrix(NA, nrow = 2, ncol = T + 1)
  V  <- array(NA, dim = c(2, 2, T + 1))
  
  # Initial state: (y1, y0) with independent normals
  mu[,1] <- c(m0, m0)
  V[,,1] <- diag(c(s0^2, s0^2))
  
  # Recurse forward
  for (t in 2:(T + 1)) {
    mu[,t] <- A %*% mu[,t-1] + c
    V[,,t] <- A %*% V[,,t-1] %*% t(A) + Sigma_u
  }
  
  # Extract marginal means and variances of y_t
  y_means <- mu[1,]
  y_vars  <- V[1,1,]
  
  list(
    A = A,
    y_means = y_means,
    y_sds   = sqrt(y_vars),
    stationary = all(Mod(eigen(A)$values) < 1)
  )
}

# ---------- UI ----------
ui <- fluidPage(
  
  titlePanel("Gaussian AR(2) marginals and paths"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("b0",
                  "β₀",
                  min = -5,
                  max = 5,
                  value = 0,
                  step = 0.1),
      sliderInput("b1",
                  "β₁",
                  min = -2,
                  max = 2,
                  value = 0,
                  step = 0.05),
      sliderInput("b2",
                  "β₂",
                  min = -2,
                  max = 2,
                  value = 0,
                  step = 0.05),
      sliderInput("s",
                  "σ",
                  min = 0,
                  max = 2,
                  value = 1, 
                  step = 0.1),
      sliderInput("m0",
                  "μ₀",
                  min = -5,
                  max = 5,
                  value = 0,
                  step = 0.1),
      sliderInput("T",
                  "T",
                  min = 20,
                  max = 200,
                  step = 20,
                  value = 100),
      actionButton("redo", "New sample path"),
    ),
    
    mainPanel(
      plotOutput("fanChart", height = "400px"),
      plotOutput("stationarityPlot", height = "400px"),
      verbatimTextOutput("stationarity")
    )
  )
)

# ---------- Server ----------
server <- function(input, output) {
  
  output$fanChart <- renderPlot({
    input$redo
    
    b0 <- input$b0
    b1 <- input$b1
    b2 <- input$b2
    s  <- input$s
    m0 <- input$m0
    s0 <- 1
    T  <- input$T
    
    # Compute state-space means and variances
    out <- ar2_state_space(T, b0, b1, b2, s, m0, s0)
    
    range <- 0:T
    alpha <- c(0.01, seq(0.1, 0.9, by = 0.1))
    
    middle <- out$y_means
    sds    <- out$y_sds
    
    plot(range, middle, type = "l",
         xaxt = "n", 
         yaxt = "n",
         xlab = "t",
         ylab = expression(y[t]),
         ylim = c(-20, 20), bty = "n",
         col = "white")
    
    for(a in alpha){
      U <- qnorm(1 - a / 2, mean = middle, sd = sds)
      L <- qnorm(a / 2, mean = middle, sd = sds)
      
      polygon(
        c(range, rev(range)),
        c(U, rev(L)),
        col = rgb(1, 0, 0, 0.15),
        border = NA
      )
    }
    
    inc <- 20
    axis(1, pos = 0, at = seq(0, max(range), by = inc), 
         labels = c(NA, seq(inc, max(range), by = inc)))
    axis(2, pos = 0)
    
    lines(range, simulate_ar_2(max(range) + 1, b0, b1, b2, s, m0, s0), 
          col = "black", lwd = 2)
  })
  
  output$stationarityPlot <- renderPlot({
    b1 <- input$b1
    b2 <- input$b2
    
    # triangle vertices
    tri_x <- c(2, -2, 0)
    tri_y <- c(-1, -1, 1)
    
    plot(tri_x, tri_y, type = "n",
         xlab = expression(beta[1]),
         ylab = expression(beta[2]),
         xlim = c(-2.5, 2.5),
         ylim = c(-2.5, 1.5),
         xaxt = "n", yaxt = "n", bty = "n")
    axis(1, at = c(-2, 0, 2), pos = 0)
    axis(2, at = c(-1, 0, 1), pos = 0)
    
    polygon(tri_x, tri_y, col = rgb(0, 0, 1, 0.2), border = "blue", lwd = 2)
    points(b1, b2, pch = 19, col = "red", cex = 1.5)
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
