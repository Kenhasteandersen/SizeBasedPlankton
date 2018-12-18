#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
require(latex2exp)
source("basetools.R")
source("model.R")

cex = 1 #  Size of text labels on plots


# ===================================================
# Define UI for application
# ===================================================

ui <- fluidPage(
  
  # Application title
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
                  "Light (PAR; mu mol photons/m2/s)",
                  min = 0,
                  max = 500,
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
                  "Diffusion rate (m/day)",
                  min = 0,
                  max = 5,
                  step = 0.1,
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
                    width="700px")),
      plotOutput("plotSpectrum", width="700px"),
      plotOutput("plotRates", width="700px"),
      plotOutput("plotFunction", width="600px"),
      plotOutput("plotTime", width="700px")#,
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
    # Simulate:
    nSave = 1000
    
    tEnd = p$tEnd
    if (p$amplitudeL>0)
      tEnd = 2*365
    
    out = ode(c(0.1*p$N0, p$DOC0, p$POM0, p$B0), seq(0, tEnd, length.out = nSave), derivative, p)
    nSave = dim(out)[1]
    # Assemble results:
    ix = seq(floor(nSave/2),nSave)
    ixB = 5:(p$n+4)
    
    Bmin = 0*p$m
    Bmax = 0*p$m
    for (i in 1:p$n) {
      Bmin[i] = min(out[ix,ixB[i]])
      Bmax[i] = max(out[ix,ixB[i]])
    }
    
    result = list(
      p = p,
      t = out[,"time"],
      y = out[,2:(p$n+4)],
      
      N = mean(out[ix,2]),
      DOC = mean(out[ix,3]),
      POM = mean(out[ix,4]),
      B = colMeans(out[ix,ixB]),
      
      Bmin = Bmin,
      Bmax = Bmax)
    
    result = c(result, list(rates = calcRates(max(result$t), result$N, result$DOC, result$POM, result$B,p)))
    return(result)
    
  })
  #
  # Plot spectrum
  #
  output$plotSpectrum <- renderPlot({
    p = sim()$p
    m = p$m
    
    par(cex.axis=cex,
        cex.lab=cex,
        mar=c(4, 5, 6, 2) + 0.1)
    
    alpha = 0.25
    colOsmo = rgb(0.5,0,0.5,alpha=alpha)
    colPhoto = rgb(0,1,0,alpha=alpha)
    colN = rgb(0,0,1,alpha=alpha)
    colMixo = rgb(1,0.5,0.5,alpha=alpha)
    colHetero = rgb(1,0,0,alpha=alpha)
    
    ylim = c(5,1000)
    fac = sqrt(m[2]/m[1])
    
    ixt = which(floor(sim()$t)==input$t+365)[1]
    if (input$amplitudeL==0) {
      B = sim()$B
      N = sim()$N
      DOC = sim()$DOC
      r = sim()$rates
    } else {
      B = sim()$y[ixt, 4:(p$n+3)]
      N = sim()$y[ixt, 1]
      DOC = sim()$y[ixt,2]
      r = calcRates(input$t, N, DOC, 0, B, p)
    }
    
    plot(m, B, 
         log="xy", type="b", lwd=8,
         ylim=ylim, 
         xlab=TeX("Carbon mass ($\\mu$gC)"),
         ylab=TeX("Biomass ($\\mu$gC/l)"),
         mar=c(4,5,8,2)+0.1)
    if (input$amplitudeL==0)
      polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim()$Bmin, sim()$Bmax[seq(p$n,1,by = -1)]), 
              col=rgb(0.5,0.5,0.5,alpha=alpha), border=NA)
    
    # Determine limiting process:
    ixOsmotroph = (r$JDOCreal > r$JLreal)
    ixLlimited = ((r$JCtot < r$JNtot*p$rhoCN) & !ixOsmotroph)
    ixHetero = (r$JNloss_feeding>0)
    ixMixo = ((r$JFreal/r$JLreal > 0.25) & !ixHetero)
    ixMixo[is.na(ixMixo)] = 0
    
    for (i in 1:p$n) {
      col = colN
      if (ixOsmotroph[i])
        col = colOsmo  
      if (ixHetero[i])
        col = colHetero
      if (ixMixo[i])
        col = colMixo
      if ((!ixMixo[i]) & !ixHetero[i] & ixLlimited[i])
        col = colPhoto
      
      #if ((!ixPhototroph[i]) & ixMixo[i])
      #  col = colMixo
      polygon(m[i]*c(1/fac,fac,fac,1/fac), c(0.5*ylim[1], 0.5*ylim[1], 2*ylim[2], 2*ylim[2]),
              col=col,
              lty=0)
    }
    
    # Plot biomass
    #lines(m,sim()$B, lwd=8, type="b")
    # Plot nutrients
    #lines(c(0.5*m[1], 2*m[p$n]), 1000*sim()$N*c(1,1), lwd=3, col="green")
    
    #
    # Add extra size labels
    #
    d = 10^seq(-6,-1,by=1)
    axis(side=3, line=0,
         at=0.3e6*d^3,
         labels = c("0.01","0.1","1","10","100","1e3"))
    mtext(TeX("Diameter ($\\mu$m)"), side=3, line=3, at=1e-3, adj=1,cex=cex)
    
    legend(x="topright", bty="n", cex=cex,
           legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
           fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
           border=c("black","black","black","black","black","transparent"),
           lwd = c(0,0,0,0,0,3))
    
    text(x=m[1], y=10, labels=TeX(sprintf("DIN: %2.2f $\\mu$gN/l", N)) , cex=cex, pos=4, col=grey(0.5))
    text(x=m[1], y=7, labels=TeX(sprintf("DOC: %2.2f $\\mu$gC/l", DOC)), cex=cex, pos=4, col=grey(0.5))
    
    func = calcFunctions(sim()$p, sim()$rates, sim()$N, sim()$B)
    text(x=1e-2, 14, 
         labels=TeX(sprintf("Picoplankton: %2.2f $gC/m$^2$", func$Bpico)),
         cex=cex, pos=4, col=grey(0.5))
    text(x=1e-2, 10, 
         labels=TeX(sprintf("Nanoplankton: %2.2f $gC/m$^2$", func$Bnano)),
         cex=cex, pos=4, col=grey(0.5))
    text(x=1e-2, 7, 
         labels=TeX(sprintf("Microplankton: %2.2f $gC/m$^2$", func$Bmicro)),
         cex=cex, pos=4, col=grey(0.5))
    #text(x=m[1], y=7 , labels=sprintf("POM: %2.2f mugC/l",sim()$POM ), cex=cex, pos=4, col=grey(0.5))
    
    box()
    
  })
  #
  # Plot rates:
  #
  output$plotRates <- renderPlot({
    p = sim()$p
    
    mm = 10^seq(-8,2, length.out = 100)  
    
    ixt = which(floor(sim()$t)==input$t+365)[1]
    if (input$amplitudeL==0) {
      L = input$L
      B = sim()$B
      N = sim()$N
      DOC = sim()$DOC
      r = sim()$rates
    } else {
      L = input$L*0.5*(1+input$amplitudeL*(-cos(input$t*2*pi/365))) 
      B = sim()$y[ixt, 4:(p$n+3)]
      N = sim()$y[ixt, 1]
      DOC = sim()$y[ixt,2]
      r = calcRates(input$t, N, DOC, 0, B, p)
    }
    
    par(cex.axis=cex,
        cex.lab=cex,
        mar=c(4, 5, 6, 2) + 0.1)
    ylim = c(0,1.5)
    # plot(mm, p$AL*mm^(2/3)*input$L/mm, lwd=4, type="l", col="green", log="x", xlim=range(p$m),
    #      ylim=ylim,
    #      xlab="Carbon mass (mu gC)",
    #      ylab="Rates (1/day)")
    # 
    # lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim()$N/mm, lwd=4, col="blue")
    # lines(mm, p$AN*mm^(1/3)*sim()$DOC/mm, lwd=4, col="brown")
    # #lines(mm, AF*mm*mean(r$F)/mm, lty=3, lwd=1, col="red")
    # lines(p$m, p$AFm*r$F/p$m, lwd=4, col="red")
    plot(p$m, p$Jmax/p$m*r$f, lwd=10, type="l", col="black", log="x", xlim=range(p$m),
         ylim=ylim,
         xlab=TeX("Carbon mass ($\\mu$gC)"),
         ylab="Rates (1/day)")
    lines(p$m, p$Jmax/p$m, lty=3)
    
    #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
    lines(p$m, p$ALm*input$L/p$m, lty=3, lwd=1, col="green")
    lines(p$m, r$JLreal/p$m, lwd=4, col="green")
    
    lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim()$N/mm, lwd=1, lty=3, col="blue")
    lines(p$m, r$JNreal/p$m*p$rhoCN, lwd=4, col="blue")
    
    lines(mm, p$AN*mm^(1/3)*sim()$DOC/mm, lwd=1, lty=3, col="brown")
    lines(p$m, r$JDOCreal/p$m, lwd=4, col="brown")
    
    lines(p$m, r$JFreal/p$m,lwd=4,col="red")
    
    legend(x="topright", cex=cex,
           legend=c("Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
           col=c("green","blue","brown","red","black"),
           lwd=4,
           bty="n")
    # plot(m, r$dBdt/B, log="x", type="l", lwd=2, ylim=c(-0.2,1))
    #  lines(m, r$mortpred, col="red")
    #  lines(m, r$mort, col="red", lty=2)
    #  lines(m,0*m,lty=3)
  })
  #
  # Complex rates plot
  #
  output$plotComplexRates <- renderPlot({
    p = sim()$p
    
    ixt = which(floor(sim()$t)==input$t+365)[1]
    if (input$amplitudeL==0) {
      L = input$L
      B = sim()$B
      N = sim()$N
      DOC = sim()$DOC
      r = sim()$rates
    } else {
      L = input$L*0.5*(1+input$amplitudeL*(-cos(input$t*2*pi/365))) 
      B = sim()$y[ixt, 4:(p$n+3)]
      N = sim()$y[ixt, 1]
      DOC = sim()$y[ixt,2]
      r = calcRates(input$t, N, DOC, 0, B, p)
    }
    
    par(cex.axis=cex,
        cex.lab=cex,
        mar=c(4, 5, 6, 2) + 0.1)

    m = p$m
    plot(m, p$rhoCN*r$JN/m, log="x", type="l", col="blue", 
         ylim=c(-1,2.5), xlim=c(1e-8, 1000),
         lwd=3,
         xlab="Carbon mass (mu gC)",
         ylab="Rates (1/day)")
    lines(m, r$JL/m, col="green", lwd=3)
    lines(m, r$JF/m, col="red", lwd=3)
    lines(m, r$Jtot/m, lty=2)
    lines(m, p$Jmax/m*r$f)
    lines(m, r$dBdt/B,lwd=2)
    lines(m,0*m,lty=3)
    lines(m, -r$mortpred, col="red")
    lines(m, -r$mort, col="red", lty=2)
    lines(m, -p$Jresp/m, col="magenta")
    lines(m, -p$mort2*B, col="blue")
    lines(m, -p$mortHTL*(m>=p$mHTL), col="orange")
    legend(x="topright",
           legend=c("rhoCN*JN/m","JL/m","JF/m","Jtot/m","Jmax/m","r","0","mort_pred","mort","resp","mort2","mortHTL"),
           col=c("blue", "green", "red", "black","black","black","black","red", "red", "magenta", "blue","orange"),
           lty=c(1,1,1,2,1,1,3,1,2,1,1,1),
           lwd=c(3,3,3,1,1,2,1,1,1,1,1,1))
    
    
    
  })
  
  #
  # Plot time-line
  #
  output$plotTime <- renderPlot({
    p = sim()$p
    t = sim()$t
    
    par(cex.axis=cex,
        cex.lab=cex,
        mar=c(4, 5, 6, 2) + 0.1)
    
    y = sim()$y
    y[y <= 0] = 1e-30
    
    ylim = c(max(1e-5, min(sim()$y)), max(sim()$y))
    if (input$amplitudeL==0) {
      xlim = range(t)  
    } else {
      xlim = c(0,365)
      t = 365+t-max(t)
    }
    
    plot(t, y[,2], log="y", type="l", col="magenta", 
         ylim=ylim, xlim=xlim, lwd=2,
         xlab="Time (day)", ylab=TeX("Biomass ($\\mu$gC/l)"))
    lines(t, y[,1], col="blue", lwd=2)
    lines(t, y[,3], col="orange", lwd=2)
    for (i in 1:p$n)
      lines(t, y[,i+3], lwd=i/p$n*3, col="black")
    
    if (input$amplitudeL>0) {
      lines((input$t)*c(1,1), ylim, lty=3)
      lines(t, input$L*0.5*(1+input$amplitudeL*(-cos(t*2*pi/365))),
            col="orange", lwd=2)
    }
  })
  #
  # Plot functions:
  #
  output$plotFunction <- renderPlot({
    func = calcFunctions(sim()$p, sim()$rates, sim()$N, sim()$B)
    par(mar=c(5,12,4,2))
    barplot(height=c(func$prodNew, func$prodCgross, func$prodCnet, func$prodHTL),
            names.arg = c("New production", "Gross PP", "Net PP", "HTL"),
            xlab = TeX("Production (gC/m$^2$/yr)"),
            horiz=TRUE, las=1,
            border=NA)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

