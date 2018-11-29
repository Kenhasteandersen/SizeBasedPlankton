#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  h1('Size-based plankton simulator'),
  p('Simulate a plankton ecosystem in the upper part of a watercolumn. 
   Cell size is the only trait characterizing each plankton group.
    All groups are able to perform photoharvesting, taking up inorganic nutrients, and do phagotrophy.
    The trophic strategy is an emergent property.')
  ,
  # Sidebar with a slider input for number of bins 
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
                  max = 5,
                  step = 0.1,
                  value = 1)
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
                  max = 0.01,
                  step = 0.0002,
                  value = 0.001)
      #,
      #sliderInput("AF",
      #             "AF (mugC^(1/3)/(W/m^2)^(-1)/d",
      #             min = 0,
      #             max = 0.02,
      #             value = 0.006,
      #             step = 0.001)
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("Plots",
                 height="1000px")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  output$Plots <- renderPlot({
    phi = function(z, beta=500, sigma=1.3)
      exp( -(log(z/beta))^2/(2*sigma^2) )
    
    calcrates = function(t,N,DOC,POM,B) {
      JN =   AN*m^(1/3)*N  # Diffusive nutrient uptake
      JDOC = AN*m^(1/3)*DOC # Diffusive DOC uptake
      JL =   AL*m^(2/3)*L  # Photoharvesting
      
      for (i in 1:n) # loop over predators
        F[i] = sum(theta[i,]*B)
      JF = AF*m*F        # Feeding xw
      
      JCtot = JL+epsilon*JF-Jresp+JDOC # Total carbon intake
      JNtot = JN*rhoCN+epsilon*JF # In units of carbon
      Jtot = pmin( JCtot, JNtot )  # Liebig law
      
      f = (Jtot) / (Jtot + Jmax) # feeding level
      JNloss = (1-f)*(JNtot-Jtot)
      JCloss = (1-f)*(JCtot-Jtot)
      
      mortpred = rep(0,n)
      for (j in 1:n)  # loop over prey
        mortpred[j] = sum(theta[,j]*(1-f)*AF*B)
      
      dBdt = Jmax*f*B/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL))*B
      dNdt = d*(N0-N) - sum(f*JN*B/m) + sum(JNloss*B/m)
      dDOCdt = sum(JCloss*B/m) -sum(f*JDOC*B/m) - d*DOC
      dPOMdt = 0#(1-epsilon)*sum(f*JF*B/m) - d*POM
      
      return(list(m=m, 
                  N=N, B=B, DOC=DOC, POM=POM,
                  dNdt=dNdt, dDOCdt=dDOCdt, dPOMdt=dPOMdt, dBdt=dBdt, 
                  JN=JN, JL=JL, JF=JF, Jtot=Jtot, f=f, F=F,
                  mortpred=mortpred, mort=mort,
                  totKilled = sum(JF*(1-f)*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
    }
    
    plotrates = function(t,N,DOC,POM,B) {
      r = calcrates(t,N,DOC,POM,B)
      m = r$m
      
      par(mfcol=c(3,1),
          cex.axis=2,
          cex.lab=2,
          mar=c(4, 5, 6, 2) + 0.1)
      
      alpha = 0.25
      colPhoto = rgb(0,1,0,alpha=alpha)
      colN = rgb(0,0,1,alpha=alpha)
      colMixo = rgb(1,0.5,0.5,alpha=alpha)
      colHetero = rgb(1,0,0,alpha=alpha)
      
      ylim = c(5,1000)
      fac = sqrt(m[2]/m[1])
      plot(m, r$B, 
           log="xy", type="b", lwd=4,
           ylim=ylim, 
           xlab="Carbon mass (mu gC)",
           ylab="Biomass (mu gC/l)",
           mar=c(4,5,8,2)+0.1)
      # Determine limiting process:
      ixPhototroph = (r$Jtot > 5*r$JF)
      ixNlimited = (r$JL>rhoCN*r$JN)
      ixMixo = (r$JL > 0.2*r$JF)
      for (i in 1:n) {
        col = colHetero
        if (ixPhototroph[i]) 
          col = colPhoto
        if (ixPhototroph[i] & ixNlimited[i])
          col = colN
        if ((!ixPhototroph[i]) & ixMixo[i])
          col = colMixo
        polygon(m[i]*c(1/fac,fac,fac,1/fac), c(0.5*ylim[1], 0.5*ylim[1], 2*ylim[2], 2*ylim[2]),
                col=col,
                lty=0)
      }
      text(x=m[1], y=10, labels=sprintf("DIN: %2.2f mugN/l",N ), cex=2, pos=4)
      text(x=m[1], y=6 , labels=sprintf("DOC: %2.2f mugC/l",DOC ), cex=2, pos=4)

      # Plot biomass
      lines(m,r$B, lwd=8, type="b")
      # Plot nutrients
      lines(c(0.5*m[1], 2*m[n]), 1000*r$N*c(1,1), lwd=3, col="green")
      
      legend(x="topright", bty="n", cex=2,
             legend=c("Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs","N (pgN/l)"),
             fill=c(colPhoto,colN,colMixo,colHetero,"transparent"),
             border=c("black","black","black","black","transparent"),
             lwd = c(0,0,0,0,3),
             col = c("green","blue","pink","red","green"))
             #border=c("black","black","black","black","red"))
      box()
      #
      # Add extra size labels
      #
      l = 10^seq(-6,-1,by=1)
      axis(side=3, line=0,
           at=0.3e6*l^3,
           labels = c("1e-6","1e-5","1e-4","1e-3","1e-2","1e-1"))
      mtext("Diameter (cm)", side=3, line=3, at=1e-3, adj=1,cex=1.25)
      
      #
      # Simple rates-plot
      #
      mm = 10^seq(-8,2, length.out = 100)  
      ylim = c(0,1.5)
      plot(mm, AL*mm^(2/3)*input$L/mm, lwd=4, type="l", col="green", log="x", xlim=range(m),
           ylim=ylim,
           xlab="Carbon mass (mu gC)",
           ylab="Rates (1/day)")
      
      lines(mm, AN*mm^(1/3)*rhoCN*r$N/mm, lwd=4, col="blue")
      lines(mm, AN*mm^(1/3)*r$DOC/mm, lwd=4, col="brown")
      #lines(mm, AF*mm*mean(r$F)/mm, lty=3, lwd=1, col="red")
      lines(m, AF*m*r$F/m, lwd=4, col="red")
      lines(m, Jmax/m*r$f, lwd=8)
      legend(x="topright", cex=2,
             legend=c("Light uptake","Nutrient uptake","DOM uptake","Food consumption","Division rate"),
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
    }
    
    derivative = function(t,y,par) {
      N = y[1]
      DOC = y[2]
      POM = y[3]
      B = y[4:(3+n)]
      
      rates = calcrates(t,N,DOC,POM,B)
      
      return(list(c(rates$dNdt, rates$dDOCdt, rates$dPOMdt, rates$dBdt)))
    }
    
    n = 9; # No of groups
    m = 10^seq(-7,1,length.out = n)  # mu g N
    mHTL = m[n-2]
    
    rhoCN = 5.68 # C:N ratio
    epsilon = 0.7 # Assimilation efficiency
    #
    # Clearance rates:
    #
    #factor = 1.5*(1e-6*rhoCN)^(1/3) 
    AN = 0.000162# Again, mathilde.  2  2.5e-3*factor# l/day/mugN^1/3
    AL = 0.01 # Following Mathildes defense paper # 72*factor^2#  mugC/day/(W m^2)/mugC^(2/3)
    AF = 0.006 #input$AF #0.0124 # Again, Mathilde # 1.3e3*factor^3 # l/d/mugN
    
    #
    # Prey encounter
    #
    theta = matrix(nrow=n, ncol=n)
    for (i in 1:n)
      theta[i,] = phi(m[i]/m)   # (predator, prey)
    
    Jmax = 1 * m # mugC/day
    Jresp = 0.04*Jmax
    mort = 0*0.005*(Jmax/m) * m^(-1/4)
    mort2 = input$mort2 #0.0015
    mortHTL = input$mortHTL
    
    d = input$d # diffusion rate
    N0 = 150#input$N0 # Deep nutrient levels
    DOC0 = 0
    POM0 = 0
    B0 = rep(10,n)
    L = input$L
    
    
    
    # plot(m, rhoCN*AN*m^(1/3)*4, log="xy", type="b", xlab="mugN", ylab="mugC/day", ylim=c(1e-7,1e-1))
    # lines(m, AL*m^(2/3)*2, col="green")
    # lines(m, AF*m*5, col="red")
    
    # l = 1.5*(1e-6*m*rhoCN)^(1/3)
    # plot(l, 72*l^2*2, log="xy", type="l", col="green", ylim=c(1e-7,1e-1))
    # lines(l, rhoCN*2.5e-3*l^1*4,col="blue")
    # lines(l, 1.3e3*l^3*5,col="red")
    
    nSave = 100
    out = ode(c(N0,DOC0, POM0, B0), seq(0, 300, length.out = nSave), derivative)
    ixB = 5:(n+4)
    
    #par(mfcol=c(2,1))
    #t = out[,"time"]
    #plot(t, out[,"1"], log="y", ylim=c(1e-5,1000), type="l", col="blue")
    #for (j in 1:n)
    #  lines(t, out[,j+2], lwd=j/3)
    
    #plot(m, out[nSave, 3:(n+2)], type="b", log="xy")
    
    #r=calcrates(t[nSave], out[nSave,2], out[nSave,ixB])
    plotrates(t[nSave], out[nSave,2], out[nSave,3], out[nSave,4], out[nSave,ixB])
  })
}

library(deSolve)



# Run the application 
shinyApp(ui = ui, server = server)

