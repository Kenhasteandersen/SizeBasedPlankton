#
# Install app:
#  ssh ken@oceanlife.dtuaqua.dk
#  update the git (git pull)
#  cp ~/SizeBasedPlankton/R/* to /srv/shiny-server/Watercolumn
#  cp ~/SizeBasedPlankton/Cpp/model.cpp to /srv/shiny-server/Watercolumn
#  If using new packages install them by running R as root (sudo su; R; install.packages("XXX))
#  sudo systemctl restart shiny-server
#

library(shiny)
options(shiny.sanitize.errors = FALSE)

source("modelWatercolumn.R")

# ===================================================
# Define UI for application
# ===================================================

uiWatercolumn <- fluidPage(
  tags$head(
    # Add google analytics tracking:
    includeHTML(("googleanalyticsWatercolumn.html")),
    # Make rules widers:
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  )
  ,
  h1('Size-based water column plankton simulator'),
  p('Simulate a plankton ecosystem in a watercolumn. 
    Cell size is the only trait characterizing each plankton group.
    All groups are able to perform photoharvesting, taking up dissolve nutrients and carbon, and do phagotrophy.
    The trophic strategy is an emergent property.'),
  HTML('<p>Made in R with shiny. 
    Code on <a href="https://github.com/Kenhasteandersen/SizeBasedPlankton">github</a>.</p>'),
  HTML('<p>THIS IS WORK IN PROGRESS. Version 0.1. January 2020. 
    See also <a href="http://oceanlife.dtuaqua.dk/Plankton">chemostat version</a>
       and other <a href="http://oceanlife.dtuaqua.dk">apps</a>. 
       By <a href="mailto:kha@aqua.dtu.dk">Ken H Andersen</a>.</p>')
  ,
  # Sidebar with a slider inputs
  sidebarLayout(
    sidebarPanel(
      HTML("<h3>Environment:</h3>"),
      sliderInput("Lsurface",
                  "Surface light (PAR; uE/m2/s)",
                  min = 0,
                  max = 300,
                  step=1,
                  value = parametersWatercolumn()$Lsurface)
      ,
      sliderInput("diff",
                  "Diffusion (m2/day)",
                  min = 0,
                  max = 50,
                  step = 1,
                  value = 1)
      ,
      sliderInput("T",
                  "Temperature",
                  min = 0,
                  max = 25,
                  step = 0.5,
                  value = 10)
      ,
      sliderInput("N0",
                  "Bottom nutrient concentration",
                  min = 0,
                  max = 200,
                  step=1,
                  value = 150)
      ,
      sliderInput("Depth",
                  "Depth (m)",
                  min=10,
                  max=200,
                  step=1,
                  value=parametersWatercolumn()$depth)
      ,
      hr(),
      HTML("<h3>Parameters:</h3>")
      ,
      sliderInput("mortHTL",
                  "Mortality by higher trophic levels on the largest sizes (1/day)",
                  min = 0,
                  max = 0.5,
                  step = 0.01,
                  value = 0.05)
      ,
      sliderInput("mort2",
                  "Quadratic mortality coef. (virulysis; 1/day/mugC)",
                  min = 0,
                  max = 0.001,
                  step = 0.00001,
                  value = 0.0002)
      ,
      sliderInput("epsilon_r",
                  "Remineralisation of virulysis and HTL losses",
                  min = 0,
                  max = 1,
                  value = 0,
                  step = 0.1)
      ,
      hr(),
      HTML("<h3>Numerical parameters</h3>")
      ,
      sliderInput("dt",
                  "Time step (days)",
                  min = 0.02,
                  max = 0.5,
                  step=0.01,
                  value=0.02)
      ,
      sliderInput("nGrid",
                  "No. of grid points",
                  min=10,
                  max=200,
                  step=1,
                  value=50)
    ),
    #
    # Show plots:
    #
    mainPanel(
      #
      # Main output panel
      #
      tabsetPanel(
        tabPanel(
          "Main output",
          plotOutput("plotWatercolumn", width="600px", height="400px"),
          plotOutput("plotFunction", width="600px", height="300px")
        )
        ,
        tabPanel(
          "Timeline",
          plotOutput("plotWatercolumnTime", width="600px")
        )
      )
    )))
#
# Define server logic
#
serverWatercolumn <- function(input, output) {
  #
  # Simulate the system when a parameter is changed:
  #
  sim <- eventReactive({
    input$Lsurface
    input$diff
    input$T
    input$N0
    input$Depth
    input$mort2
    input$mortHTL
    input$epsilon_r
    input$dt
    input$nGrid
  },
  {
    # get all base parameters
    p <- parametersWatercolumn() 
    # Update parameters with user input:
    p$Lsurface = input$Lsurface
    p$diff = input$diff
    p$T = input$T
    p$N0 = input$N0
    p$depth = input$Depth
    p$mort2 = input$mort2*p$n
    p$mortHTL = input$mortHTL
    p$remin = input$epsilon_r
    p$dt = input$dt
    p$nGrid = input$nGrid
    # Simulate
    return(simulateWatercolumn(p))
  })
  #
  # Plots:
  #
  output$plotWatercolumn <- renderPlot(plotWatercolumn(sim()))
  output$plotWatercolumnTime <- renderPlot(plotWatercolumnTime(sim()))
  output$plotFunction <- renderPlot(plotFunctionsWatercolumn(sim()))
}

# Run the application 
shinyApp(ui = uiWatercolumn, server = serverWatercolumn)

