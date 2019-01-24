#
# Install app:
#  ssh ken@172.23.130.112:

library(shiny)
options(shiny.sanitize.errors = FALSE)

source("model.R")
source("plots.R")

# ===================================================
# Define UI for application
# ===================================================

ui <- fluidPage(
  h1('Size-based plankton simulator'),
  p('Simulate a plankton ecosystem in the upper part of a watercolumn. 
   Cell size is the only trait characterizing each plankton group.
    All groups are able to perform photoharvesting, taking up dissolve nutrients and carbon, and do phagotrophy.
    The trophic strategy is an emergent property.'),
  p('THIS IS WORK IN PROGRESS. Version 0.5. November 2018.')
  ,
  # Sidebar with a slider inputs
  sidebarLayout(
    sidebarPanel(
      sliderInput("L",
                  "Light (PAR; uE/m2/s)",
                  min = 0,
                  max = 70,
                  step=1,
                  value = parameters()$L)
      ,
      sliderInput("amplitudeL",
                  "Seasonal variation (experimental)",
                  min = 0,
                  max = 1,
                  step=0.1,
                  value = 0)
      ,
      #sliderInput("N0",
      #            "Deep nutrients (mugN/l)",
      #            min = 0,
      #            max = 200,
      #            step = 0.1,
      #            value = 14)
      sliderInput("d",
                  "Exchange rate (m/day)",
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = parameters()$d)
      ,
      sliderInput("M",
                  "Thickness of productive layer (m)",
                  min = 0,
                  max = 100,
                  step = 1,
                  value = parameters()$M)
      ,
      sliderInput("mortHTL",
                  "Mortality by higher trophic levels on the largest sizes (1/day)",
                  min = 0,
                  max = 0.5,
                  step = 0.01,
                  value = parameters()$mortHTL)
      ,
      sliderInput("mort2",
                  "Quadratic mortality coef. (1/day/mugC)",
                  min = 0,
                  max = 0.001,
                  step = 0.00001,
                  value = 0.0002)
      #,
      #sliderInput("AF",
      #             "AF (mugC^(1/3)/(W/m^2)^(-1)/d",
      #             min = 0,
      #             max = 0.02,
      #             value = 0.006,
      #             step = 0.001)
      
    ),
    #
    # Show plots:
    #
    mainPanel(
      conditionalPanel(
        condition = "input.amplitudeL>0",
        sliderInput("t",
                    "Time (day)",
                    min=0, 
                    max=365,
                    step=1,
                    value=120,
                    width="600px")),
      plotOutput("plotSpectrum", width="600px", height="300px"),
      plotOutput("plotRates", width="600px", height="300px"),
      plotOutput("plotLeaks", width="600px", height="200px"),
      plotOutput("plotFunction", width="600px"),
      plotOutput("plotTime", width="600px")#,
      #plotOutput("plotComplexRates", width="700px")
    )
  )
)
#
# Define server logic
#
server <- function(input, output) {
  
  #
  # Simulate the system when a parameter is changed:
  #
  sim <- eventReactive({
    input$L
    input$amplitudeL
    input$d
    input$M
    input$mort2
    input$mortHTL
  },
  {
    # get all base parameters
    p <- parameters() 
    # Update parameters with user input:
    p$L = input$L
    p$amplitudeL = input$amplitudeL
    p$d = input$d
    p$mort2 = input$mort2*p$n
    p$mortHTL = input$mortHTL
    p$M = input$M

    if (p$amplitudeL>0)
      p$tEnd = 2*365
    
    # Simulate
    return(simulate(p))   
  })
  #
  # Plots:
  #
  output$plotSpectrum <- renderPlot(plotSpectrum(sim(), input$t))
  output$plotRates <- renderPlot(plotRates(sim(), input$t))
  output$plotLeaks = renderPlot(plotLeaks(sim(), input$t))
  output$plotComplexRates <- renderPlot(plotComplexRates(sim(), input$t))
  output$plotTime <- renderPlot(plotTimeline(sim(), input$t))
  output$plotFunction <- renderPlot(plotFunctions(sim()))
}

# Run the application 
shinyApp(ui = ui, server = server)

