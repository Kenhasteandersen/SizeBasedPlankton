#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)

cex = 1 # Size of text labels on plots

#--------------------------------------------------
# Core logic for the model
#--------------------------------------------------

parameters <- function() {
  p = list()
  
  p$n = 10; # No of groups
  p$m = 10^seq(-8,1,length.out = p$n)  # Mass bins in mugN
  p$mHTL = p$m[p$n-2] # Bins affects by HTL mortality
  
  p$rhoCN = 5.68 # C:N ratio
  p$epsilonL = 0.9
  p$epsilonF = 0.8 # Assimilation efficiency
  #
  # Clearance rates:
  #
  #factor = 1.5*(1e-6*rhoCN)^(1/3) 
  p$AN = 0.000162# Again, mathilde.  2  2.5e-3*factor# l/day/mugN^1/3
  p$AL = 0.01 # Following Mathildes defense paper # 72*factor^2#  mugC/day/(W m^2)/mugC^(2/3)
  p$AF = 0.006 #input$AF #0.0124 # Again, Mathilde # 1.3e3*factor^3 # l/d/mugN
  
  p$ANm = p$AN*p$m^(1/3)
  p$ALm = p$AL*p$m^(2/3)
  p$AFm = p$AF*p$m
  #
  # Prey encounter
  #
  p$theta = matrix(nrow=p$n, ncol=p$n)
  for (i in 1:p$n)
    p$theta[i,] = phi(p$m[i]/p$m)   # (predator, prey)
  #
  # Metabolism:
  #
  p$Jmax = 1 * p$m # mugC/day
  p$Jresp = 0.04*p$Jmax
  #
  # Losses:
  #
  p$mort = 0*0.005*(p$Jmax/p$m) * p$m^(-1/4)
  p$mort2 = 0.00015*p$n
  p$mortHTL = 0.2
  
  p$d = 1 # diffusion rate
  p$dPOM = 10 # Loss rate of POM
  p$remin = 0.4 # rate of remineralization of POM
  p$N0 = 150 # Deep nutrient levels
  p$DOC0 = 0
  p$POM0 = 0
  p$B0 = rep(10,p$n)
  p$L = 5 
  
  p$tEnd = 300 # Simulation lenght (days)
  
  return(p)
}

phi = function(z, beta=500, sigma=1.3)
  exp( -(log(z/beta))^2/(2*sigma^2) )

calcRates = function(t,N,DOC,POM,B,p) {
  with(p, {
    #
    # Potential uptakes:
    #
    JN =   ANm*N  # Diffusive nutrient uptake
    JDOC = ANm*DOC # Diffusive DOC uptake
    JL =   epsilonL*ALm*L  # Photoharvesting
    
    for (i in 1:n) # loop over predators
      F[i] = sum(theta[i,]*B)
    JF = epsilonF*AFm*F        # Feeding
    
    JCtot = JL+JF-Jresp+JDOC # Total carbon intake
    JNtot = JN+JF/rhoCN # In units of N
    Jtot = pmin( JCtot, JNtot*rhoCN )  # Liebig law; units of C
    
    f = (Jtot) / (Jtot + Jmax) # feeding level
    #
    # Actual uptakes:
    #
    JFreal = pmax(0, JF - (Jtot-f*Jmax))
    JLreal = JL-pmin(JCtot - (JF-JFreal)-f*Jmax, JL)
    JDOCreal = pmin(JDOC, Jresp + f*Jmax - JLreal - JFreal) # The min is only needed to reduce round-off errors
    JNreal = pmax(0, (f*Jmax - JFreal))/rhoCN
    
    JCcheck = f*Jmax - JLreal - JFreal - JDOCreal + Jresp
    JNcheck = f*Jmax/rhoCN - JNreal - JFreal/rhoCN - pmin(0, (f*Jmax - JFreal))/rhoCN
    # 
    # Losses:
    #
    JNloss = (1-epsilonF)/epsilonF*JFreal - pmin(0, (f*Jmax - JFreal))/rhoCN
    JCloss = (1-epsilonF)/epsilonF*JFreal/rhoCN + (1-epsilonL)/epsilonL*JLreal
    #
    # Mortality:
    #
    mortpred = rep(0,n)
    for (j in 1:n)   # loop over prey
      mortpred[j] = sum(JFreal/epsilonF*B/m * theta[,j]/rowSums(theta*B))
    #
    # System:
    #
    dBdt = Jmax*f*B/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL))*B
    dNdt   =  d*(N0-N) - sum(JNreal*B/m)   + sum(JNloss*B/m) + remin*POM/rhoCN
    dDOCdt =           - sum(JDOCreal*B/m) + sum(JCloss*B/m) + remin*POM  # 
    dPOMdt = -dPOM*POM + sum(mort2*B*B) - remin*POM
  
    return(list( 
      dNdt=dNdt, dDOCdt=dDOCdt, dPOMdt=dPOMdt, dBdt=dBdt, 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JNreal=JNreal, JDOCreal=JDOCreal, JLreal=JLreal, JFreal=JFreal,
      JNloss=JNloss, JCloss=JCloss,
      Jtot=Jtot, f=f, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      totKilled = sum(JFreal/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
  })
}

derivative = function(t,y,p) {
  N = y[1]
  DOC = y[2]
  POM = y[3]
  B = y[4:(3+p$n)]
  
  rates = calcRates(t,N,DOC,POM,B,p)
  
  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
    browser()
  
  return(list(c(rates$dNdt, rates$dDOCdt, rates$dPOMdt, rates$dBdt)))
}

# ===================================================
# Define UI for application
# ===================================================

ui <- fluidPage(
  
  # Application title
  h1('Size-based plankton simulator'),
  p('Simulate a plankton ecosystem in the upper part of a watercolumn. 
   Cell size is the only trait characterizing each plankton group.
    All groups are able to perform photoharvesting, taking up inorganic nutrients, and do phagotrophy.
    The trophic strategy is an emergent property.')
  ,
  # Sidebar with a slider inputs
  sidebarLayout(
    sidebarPanel(
      sliderInput("L",
                  "Light (W/m2)",
                  min = 0,
                  max = 5,
                  step=0.1,
                  value = 3.5)
      ,
      #sliderInput("N0",
      #            "Deep nutrients (mugN/l)",
      #            min = 0,
      #            max = 200,
      #            step = 0.1,
      #            value = 14)
      sliderInput("d",
                  "Diffusion rate (1/day)",
                  min = 0,
                  max = 1,
                  step = 0.02,
                  value = .2)
      ,
      sliderInput("mortHTL",
                  "Mortality by higher trohic levels on the three largest sizes (1/day)",
                  min = 0,
                  max = 0.5,
                  step = 0.01,
                  value = 0.2)
      ,
      sliderInput("mort2",
                  "Quadratic mortality coef. (1/day/mugC)",
                  min = 0,
                  max = 0.001,
                  step = 0.00001,
                  value = 0.0001)
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
      plotOutput("plotSpectrum", width="700px"),
      plotOutput("plotRates", width="700px"),
      plotOutput("plotTime", width="700px")
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
    input$d
    input$mort2
    input$mortHTL
  },
  {
    # get all base parameters
    p <- parameters() 
    # Update parameters with user input:
    p$L = input$L
    p$d = input$d
    p$mort2 = input$mort2*p$n
    p$mortHTL = input$mortHTL
    # Simulate:
    nSave = 1000
    out = ode(c(0.1*p$N0, p$DOC0, p$POM0, p$B0), seq(0, p$tEnd, length.out = nSave), derivative, p)
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
    
    return(list(
      p = p,
      t = out[,"time"],
      y = out[,2:(p$n+4)],
      
      N = mean(out[ix,2]),
      DOC = mean(out[ix,3]),
      POM = mean(out[ix,4]),
      B = colMeans(out[ix,ixB]),

      Bmin = Bmin,
      Bmax = Bmax)
    )
  })
  #
  # Plot spectrum
  #
  output$plotSpectrum <- renderPlot({
    p = sim()$p
    r = calcRates(sim()$t, sim()$N, sim()$DOC, sim()$POM, sim()$B, p)
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
    plot(m, sim()$B, 
         log="xy", type="b", lwd=4,
         ylim=ylim, 
         xlab="Carbon mass (mu gC)",
         ylab="Biomass (mu gC/l)",
         mar=c(4,5,8,2)+0.1)
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim()$Bmin, sim()$Bmax[seq(p$n,1,by = -1)]), 
            col=rgb(0.5,0.5,0.5,alpha=alpha), border=NA)
    
    # Determine limiting process:
    ixOsmotroph = (r$JDOC > r$JL)
    ixPhototroph = (r$Jtot > 5*r$JF)
    ixNlimited = (r$JL>p$rhoCN*r$JN)
    ixMixo = (r$JL > 0.2*r$JF)
    ixHetero = ((r$JNloss>0) & (JFreal>0))
    
    for (i in 1:p$n) {
      col = colHetero
      if (ixOsmotroph[i])
        col = colOsmo  
      if (ixPhototroph[i] & !ixOsmotroph[i]) 
        col = colPhoto
      if (ixPhototroph[i] & ixNlimited[i] & !ixOsmotroph[i])
        col = colN
      if ((!ixPhototroph[i]) & ixMixo[i])
        col = colMixo
      polygon(m[i]*c(1/fac,fac,fac,1/fac), c(0.5*ylim[1], 0.5*ylim[1], 2*ylim[2], 2*ylim[2]),
              col=col,
              lty=0)
    }
    
    # Plot biomass
    lines(m,sim()$B, lwd=8, type="b")
    # Plot nutrients
    #lines(c(0.5*m[1], 2*m[p$n]), 1000*sim()$N*c(1,1), lwd=3, col="green")
    
    legend(x="topright", bty="n", cex=cex,
           legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
           fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
           border=c("black","black","black","black","black","transparent"),
           lwd = c(0,0,0,0,0,3))

    text(x=m[1], y=14, labels=sprintf("DIN: %2.2f mugN/l",sim()$N ), cex=cex, pos=4, col=grey(0.5))
    text(x=m[1], y=10, labels=sprintf("DOC: %2.2f mugC/l",sim()$DOC ), cex=cex, pos=4, col=grey(0.5))
    text(x=m[1], y=7 , labels=sprintf("POM: %2.2f mugC/l",sim()$POM ), cex=cex, pos=4, col=grey(0.5))

    box()
    #
    # Add extra size labels
    #
    l = 10^seq(-6,-1,by=1)
    axis(side=3, line=0,
         at=0.3e6*l^3,
         labels = c("0.01","0.1","1","10","100","1e3"))
    mtext("Diameter (micro m)", side=3, line=3, at=1e-3, adj=1,cex=cex)
    
  })
  #
  # Plot rates:
  #
  output$plotRates <- renderPlot({
    p = sim()$p
    r = calcRates(sim()$t, sim()$N, sim()$DOC, sim()$POM, sim()$B, p)

    mm = 10^seq(-8,2, length.out = 100)  
    
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
         xlab="Carbon mass (mu gC)",
         ylab="Rates (1/day)")
  
    lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
    lines(p$m, r$JLreal/p$m, lwd=4, col="green")
    
    lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim()$N/mm, lwd=1, lty=3, col="blue")
    lines(p$m, r$JNreal/p$m*p$rhoCN, lwd=4, col="blue")
    
    lines(mm, p$AN*mm^(1/3)*sim()$DOC/mm, lwd=1, lty=3, col="brown")
    lines(p$m, r$JDOCreal/p$m, lwd=4, col="brown")

    lines(p$m, r$JFreal/p$m,lwd=4,col="red")

    legend(x="topright", cex=cex,
           legend=c("Light uptake","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
           col=c("green","blue","brown","red","black"),
           lwd=4,
           bty="n")
    #
    # Complex rates plot:
    #
    if (1==0) {
      plot(m, rhoCN*r$JN/m, log="x", type="l", col="blue", ylim=c(-1,1.5), lwd=3,
           xlab="Carbon mass (mu gC)",
           ylab="Rates (1/day)")
      lines(m, r$JL/m, col="green", lwd=3)
      lines(m, r$JF/m, col="red", lwd=3)
      lines(m, r$Jtot/m, lty=2)
      lines(m, Jmax/m*r$f)
      lines(m, r$dBdt/B,lwd=2)
      lines(m,0*m,lty=3)
      lines(m, -r$mortpred, col="red")
      lines(m, -r$mort, col="red", lty=2)
      lines(m, -Jresp/m, col="magenta")
      lines(m, -mort2*B, col="blue")
      lines(m, -mortHTL*(m>=mHTL), col="orange")
      legend(x="topright",
             legend=c("rhoCN*JN/m","JL/m","JF/m","Jtot/m","Jmax/m","r","0","mort_pred","mort","resp","mort2","mortHTL"),
             col=c("blue", "green", "red", "black","black","black","black","red", "red", "magenta", "blue","orange"),
             lty=c(1,1,1,2,1,1,3,1,2,1,1,1),
             lwd=c(3,3,3,1,1,2,1,1,1,1,1,1))
      
    }
    # plot(m, r$dBdt/B, log="x", type="l", lwd=2, ylim=c(-0.2,1))
    #  lines(m, r$mortpred, col="red")
    #  lines(m, r$mort, col="red", lty=2)
    #  lines(m,0*m,lty=3)
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
    
    plot(t, sim()$y[,2], log="y", type="l", col="magenta", 
         ylim=c(max(1e-5, min(sim()$y)), max(sim()$y)), lwd=2)
    lines(t, sim()$y[,1], col="blue", lwd=2)
    lines(t, sim()$y[,3], col="orange", lwd=2)
    for (i in 1:p$n)
      lines(t, sim()$y[,i+3], lwd=i/p$n*3, col="black")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

