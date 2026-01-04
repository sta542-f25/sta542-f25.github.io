library(shiny)
library(KFAS)

ui <- fluidPage(
  titlePanel("Kalman filtering and smoothing"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("true_var", "True Variance (of simulated data):",
                  min = 0.01, max = 3, value = 1, step = 0.01),
      sliderInput("error_var", "Measurement error variance:",
                  min = 0.01, max = 3, value = 1, step = 0.01),
      sliderInput("state_var", "Transition error Variance:",
                  min = 0.01, max = 3, value = 1, step = 0.01),
      checkboxInput("show_filtered", "Show Filtered Means (red)", value = TRUE),
      checkboxInput("show_smoothed", "Show Smoothed Means (blue)", value = TRUE)
    ),
    
    mainPanel(
      plotOutput("kf_plot", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  output$kf_plot <- renderPlot({
    set.seed(123)
    T <- 100
    
    # Simulate data with chosen true variance
    y <- sin((1:T) / 4) + rnorm(T, sd = sqrt(input$true_var))
    
    # Build model matrices
    Zt <- matrix(1)
    Tt <- matrix(1)
    Rt <- matrix(1)
    Qt <- matrix(input$state_var)
    Ht <- matrix(input$error_var)
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
    
    # Plot results
    plot(1:T, y, type = "l", lwd = 2, xlab = "Time", ylab = "Value",
         main = "Kalman Filter and Smoother", col = "black")
    
    if (input$show_filtered) {
      lines(1:T, c(ffbs$att), col = "red", lwd = 2)
    }
    if (input$show_smoothed) {
      lines(1:T, c(ffbs$alphahat), col = "blue", lwd = 2)
    }
    
    legend("topright", legend = c("Observed", 
                                  if (input$show_filtered) "Filtered" else NULL,
                                  if (input$show_smoothed) "Smoothed" else NULL),
           col = c("black", 
                   if (input$show_filtered) "red" else NULL,
                   if (input$show_smoothed) "blue" else NULL),
           lty = 1, lwd = 2, bty = "n")
  })
}

shinyApp(ui, server)
