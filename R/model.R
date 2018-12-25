library(deSolve)
require(latex2exp)
source("basetools.R")
library(tictoc)

#--------------------------------------------------
# Core logic for the model
#--------------------------------------------------

parameters <- function() {
  p = list()
  
  p$n = 10; # No of groups
  p$m = 10^seq(-8,1,length.out = p$n)  # Mass bins in mugC
  
  p$rhoCN = 5.68 # C:N mass ratio
  p$epsilonL = 0.9 # Light uptake efficiency
  p$epsilonF = 0.8 # Assimilation efficiency
  #
  # Cell wall fraction of mass:
  #
  p$c = 0.0015 # the constant is increased a bit to limit the lower cell size
  nu = p$c * p$m^(-1/3)  
  #
  # Clearance rates:
  #
  #factor = (1e-6)^(1/3)/1.5 
  p$AN = 0.00004 # 0.000162 # Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  #p$AL = 0.0019 # if using Al propto m^(2/3) for non-diatoms
  #p$AL = 0.0012 # if using my shading formula for non-diatoms
  #p$cL = 0.021 # if using my shading formula for non-diatoms
  p$AL = 0.000914 # if using Andys shading formula for non-diatoms
  p$cL = 21 # if using Andys shading formula for non-diatoms
  p$AF = 0.018  #  Fits to TK data for protists
  
  p$ANm = p$AN*p$m^(1/3)
  #p$ALm = p$AL*p$m^(2/3)*(1-nu)
  #p$ALm = p$cL*p$m * p$AL*p$m^(2/3) / ( p$cL*p$m + p$AL*p$m^(2/3) )  # shading formula
  p$ALm = p$AL*p$m^(2/3) * exp(1 -p$cL*p$m^(1/3) )  # shading formula
  p$AFm = p$AF*p$m
  #
  # Prey encounter
  #
  p$beta = 500
  p$sigma = 1.3
  p$theta = matrix(nrow=p$n, ncol=p$n)
  for (i in 1:p$n)
    p$theta[i,] = phi(p$m[i]/p$m, p$beta, p$sigma)   # (predator, prey)
  #
  # Metabolism:
  #
  p$alphaJ = 1
  p$Jmax = p$alphaJ * p$m * (1-nu) # mugC/day
  p$Jresp = 0.1*p$Jmax
  #
  # Losses:
  #
  p$mort = 0*0.005*(p$Jmax/p$m) * p$m^(-1/4)
  p$mort2 = 0.0002*p$n
  p$mortHTL = 0.2
  p$mHTL = max(p$m)/p$beta # Bins affecte by HTL mortality
  
  p$dPOM = 10 # Loss rate of POM
  p$epsilonPOM = 0.5 # fraction of mortality losses reminerilized to N and DOC
  
  p$d = 2 # diffusion rate
  p$M = 30 # Thickness of the mixed layer
  p$N0 = 150 # Deep nutrient levels
  
  p$DOC0 = 0
  p$POM0 = 0
  p$B0 = rep(10,p$n)
  p$L = 125  # PAR
  p$amplitudeL = 0 # amplitude of seasonal light variation in fractions of L
  
  p$tEnd = 365 # Simulation length (days)
  
  return(p)
}

phi = function(z, beta, sigma)
  exp( -(log(z/beta))^2/(2*sigma^2) )

calcRates = function(t,N,DOC,POM,B,p) {
  with(p, {
    #
    # Potential uptakes:
    #
    JN =   ANm*N  # Diffusive nutrient uptake
    JDOC = ANm*DOC # Diffusive DOC uptake
    JL =   epsilonL*ALm*L*(1+amplitudeL*(-cos(t*2*pi/365)))  # Photoharvesting
    
    F = theta %*% B
    JF = epsilonF*AFm*F        # Feeding
    
    JCtot = JL+JF-Jresp+JDOC # Total carbon intake
    JNtot = JN+JF/rhoCN # In units of N
    Jtot = pmin( JCtot, JNtot*rhoCN )  # Liebigs law; units of C
    
    f = (Jtot) / (Jtot + Jmax) # feeding level
    #
    # Actual uptakes:
    #
    JFreal = pmax(0, JF - (Jtot-f*Jmax))
    JLreal = JL-pmin((JCtot - (JF-JFreal)-f*Jmax), JL)
    JDOCreal = pmin(JDOC, Jresp + f*Jmax - JLreal - JFreal) # The min is only needed to reduce round-off errors
    JNreal = pmax(0, (f*Jmax - JFreal))/rhoCN # In units of N
    # 
    # Losses:
    #
    JNloss_feeding = -pmin(0, (f*Jmax - JFreal))/rhoCN  # Exudation of surplus nutrients by heterotrophs
    JNloss = (1-epsilonF)/epsilonF*JFreal/rhoCN + JNloss_feeding
    JCloss = (1-epsilonF)/epsilonF*JFreal + (1-epsilonL)/epsilonL*JLreal
    
    # These two should be equal to zero:
    #JCcheck = f*Jmax - JLreal - JFreal - JDOCreal + Jresp
    #JNcheck = f*Jmax/rhoCN - JNreal - JFreal/rhoCN + JNloss_feeding
    #print(JNcheck)
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JFreal/epsilonF*B/m/F)
    #
    # System:
    #
    POMgeneration = sum(mort2*B*B) + sum(mortHTL*(m>=mHTL)*B)
    
    dBdt = (Jmax*f/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL)))*B
    dNdt   =  d/M*(N0-N) - sum(JNreal*B/m)   + sum(JNloss*B/m) + epsilonPOM*POMgeneration/rhoCN
    dDOCdt =           - sum(JDOCreal*B/m) + sum(JCloss*B/m) + epsilonPOM*POMgeneration  # 
    dPOMdt = -dPOM*POM + 0*sum(mort2*B*B) + 0*sum(mortHTL*(m>=mHTL)*B) - epsilonPOM*POM
    
    # Check of nutrient conservation; should be close to zero
    #Nin = d*(N0-N)
    #Nout = (1-remin) * ( sum(mort2*B*B) + sum(mortHTL*(m>=mHTL)*B) ) / rhoCN
    #NcheckSystem = Nin - Nout - sum(dBdt)/rhoCN - dNdt
    #print(NcheckSystem)
    
    return(list( 
      dNdt=dNdt, dDOCdt=dDOCdt, dPOMdt=dPOMdt, dBdt=dBdt, 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JNreal=JNreal, JDOCreal=JDOCreal, JLreal=JLreal, JFreal=JFreal, JNloss_feeding=JNloss_feeding,
      JNloss=JNloss, JCloss=JCloss,
      Jtot=Jtot, f=f, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      POMgeneration = POMgeneration,
      totKilled = sum(JFreal/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
  })
}

derivative = function(t,y,p) {
  N = y[1]
  DOC = y[2]
  POM = y[3]
  B = y[4:(3+p$n)]
  
  rates = calcRates(t,N,DOC,POM,B,p)
  
  # Check. Expensive to evaluate, so commented out  
  #  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
  #    browser()
  
  return(list(c(rates$dNdt, rates$dDOCdt, rates$dPOMdt, rates$dBdt)))
}

calcFunctions = function(param,r,N,B) {
  with(param, {
    conversion = 365*M/10*100*1e-6 # Convert to gC/yr/m2
    #
    # New production calculated as the flux of nitrogen into the system. Units of carbon:
    #
    prodNew = conversion * d*(N0-N)/M * rhoCN  
    #
    # Primary production (carbon fixed)
    #
    prodCgross = conversion * sum(r$JLreal*B/m)/epsilonL
    prodCnet = conversion * sum(r$JLreal*B/m)
    #
    # Loss to HTL:
    #
    prodHTL = conversion * sum( mortHTL*(m>=mHTL)*B )
    #
    # Loss to depth:
    #
    prodSeq = conversion * r$POMgeneration * (1-epsilonPOM)
    #
    # Respiration
    #
    resp = conversion * sum( r$Jresp*B/m )
    #
    # Efficiencies:
    #
    effHTL = prodHTL/prodNew
    #
    # Biomasses:
    #
    conversion = M/10*100*1e-6  # Convert to gC/m2
    d = 10000 * 1.5 * (m*1e-6)^(1/3)
    Bpico = conversion * sum( B[d < 2] )
    Bnano = conversion * sum( B[d>=2 & d <20])
    Bmicro = conversion * sum( B[d>=20])
    
    return(list(
      prodNew = prodNew,
      prodCgross = prodCgross,
      prodCnet = prodCnet,
      prodHTL = prodHTL,
      resp = resp,
      effHTL = effHTL,
      Bpico=Bpico, Bnano=Bnano, Bmicro=Bmicro
    ))
  })
}


plotSpectrum <- function(sim, t=max(sim$t)) {
  p = sim$p
  m = p$m
  
#  par(cex.axis=cex,
#      cex.lab=cex,
#      mar=c(4, 5, 6, 2) + 0.1)
  
  alpha = 0.25
  colOsmo = rgb(0.5,0,0.5,alpha=alpha)
  colPhoto = rgb(0,1,0,alpha=alpha)
  colN = rgb(0,0,1,alpha=alpha)
  colMixo = rgb(1,0.5,0.5,alpha=alpha)
  colHetero = rgb(1,0,0,alpha=alpha)
  
  ylim = c(5,1000)
  fac = sqrt(m[2]/m[1])
  
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$amplitudeL==0) {
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    B = sim$y[ixt, 4:(p$n+3)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, 0, B, p)
  }
  
  defaultplot(mar=c(2.1,2.3,2.1,0))
  loglogpanel(xlim=p$m, ylim=ylim, 
              xlab="Carbon mass ($\\mu$gC)",
              ylab="Biomass ($\\mu$gC/l)")
  
  lines(m, B, type="b", lwd=8)
  #     mar=c(4,5,8,2)+0.1)
  if (p$amplitudeL==0)
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim$Bmin, sim$Bmax[seq(p$n,1,by = -1)]), 
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
  #
  # Add extra size labels
  #
  d = 10^seq(-6,-1,by=1)
  axis(side=3, line=0,
       at=0.3e6*d^3,
       labels = c("0.01","0.1","1","10","100","1e3"))
  mtext(TeX("Diameter ($\\mu$m)"), side=3, line=1.25, at=1e-3, adj=1,cex=cex)
  
  legend(x="topright", bty="n", cex=cex,
         legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
         fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
         border=c("black","black","black","black","black","transparent"),
         lwd = c(0,0,0,0,0,3))
  
  text(x=m[1], y=10, labels=TeX(sprintf("DIN: %2.2f $\\mu$mol N/l", N/14)) , cex=cex, pos=4, col=grey(0.5))
  text(x=m[1], y=7, labels=TeX(sprintf("DOC: %2.2f $\\mu$mol C/l", DOC/12)), cex=cex, pos=4, col=grey(0.5))
  
  func = calcFunctions(sim$p, sim$rates, sim$N, sim$B)
  text(x=1e-2, 14, 
       labels=TeX(sprintf("Picoplankton: %2.2f $gC/m$^2$", func$Bpico)),
       cex=cex, pos=4, col=grey(0.5))
  text(x=1e-2, 10, 
       labels=TeX(sprintf("Nanoplankton: %2.2f $gC/m$^2$", func$Bnano)),
       cex=cex, pos=4, col=grey(0.5))
  text(x=1e-2, 7, 
       labels=TeX(sprintf("Microplankton: %2.2f $gC/m$^2$", func$Bmicro)),
       cex=cex, pos=4, col=grey(0.5))
  #text(x=m[1], y=7 , labels=sprintf("POM: %2.2f mugC/l",sim$POM ), cex=cex, pos=4, col=grey(0.5))
  
  box()
}

plotRates = function(sim, t=max(sim$t)) {
  p = sim$p
  
  mm = 10^seq(-8,2, length.out = 100)  
  
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$amplitudeL==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = p$L*0.5*(1+p$amplitudeL*(-cos(t*2*pi/365))) 
    B = sim$y[ixt, 4:(p$n+3)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, 0, B, p)
  }
  
  defaultplot()
  ylim = c(-1.5,1.5)
  semilogxpanel(xlim=p$m, ylim=ylim,
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Rates (1/day)")
  #
  # Gains
  #
  lines(p$m, p$Jmax/p$m*r$f, lwd=10, type="l", col="black")# log="x", xlim=range(p$m),
  lines(p$m, p$Jmax/p$m, lty=3)
  
  #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
  lines(p$m, p$ALm*p$L/p$m, lty=3, lwd=1, col="green")
  lines(p$m, r$JLreal/p$m, lwd=4, col="green")
  
  lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim$N/mm, lwd=1, lty=3, col="blue")
  lines(p$m, r$JNreal/p$m*p$rhoCN, lwd=4, col="blue")
  
  lines(mm, p$AN*mm^(1/3)*sim$DOC/mm, lwd=1, lty=3, col="brown")
  lines(p$m, r$JDOCreal/p$m, lwd=4, col="brown")
  
  lines(p$m, r$JFreal/p$m,lwd=4,col="red")

  legend(x="topright", cex=cex,
         legend=c("Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
         col=c("green","blue","brown","red","black"),
         lwd=4,
         bty="n")
  #
  # Losses
  #
  polygon(c(1e-9,10,10,1e-9), c(-1.5,-1.5,0,0), 
          col=rgb(1,0,0,alpha=0.25), border=NA)
  
  mortHTL = p$mortHTL*(p$m>p$mHTL)
  
  lines(p$m, -(sim$rates$mortpred + mortHTL + sim$rates$mort2 + p$mort), lwd=10)
  lines(p$m, -sim$rates$mortpred, col="red", lwd=4)
  lines(p$m[p$m>=p$mHTL], -p$mortHTL*sign(p$m[p$m>p$mHTL]), col="magenta", lwd=4)
  lines(p$m, -sim$rates$mort2, col="orange", lwd=4)
  
  BSheldon =exp(mean(log(sim$B)))
  mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / (p$m[2]/p$m[1])
  lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")

  legend(x="bottomright", cex=cex,
         legend=c("Predation", "Virulysis", "Higher trophic levels"),
         col=c("red", "orange", "magenta"),
         lwd=4, bty="n")
  
  lines(p$m, 0*p$m, col="white", lwd=4)
  lines(p$m, 0*p$m, lty=dashed, lwd=2)
}

plotComplexRates = function(sim, t=max(sim$t)) {
  p = sim$p
  
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$amplitudeL==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = p$L*0.5*(1+p$amplitudeL*(-cos(t*2*pi/365))) 
    B = sim$y[ixt, 4:(p$n+3)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, 0, B, p)
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
  lines(m, 0*m,lty=3)
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
}

plotTimeline = function(sim, t=max(sim$t)) {
  p = sim$p
  t = sim$t
  
  par(cex.axis=cex,
      cex.lab=cex,
      mar=c(4, 5, 6, 2) + 0.1)
  
  y = sim$y
  y[y <= 0] = 1e-30
  
  ylim = c(max(1e-5, min(sim$y)), max(sim$y))
  if (p$amplitudeL==0) {
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
  
  if (p$amplitudeL>0) {
    lines(t*c(1,1), ylim, lty=dotted)
    lines(t, p$L*0.5*(1+p$amplitudeL*(-cos(t*2*pi/365))),
          col="orange", lwd=2)
  }
}
#
# Plot functions:
#
plotFunctions <- function(sim) {
  func = calcFunctions(sim$p, sim$rates, sim$N, sim$B)
  par(mar=c(5,12,4,2))
  barplot(height=c(func$prodNew, func$prodCgross, func$prodCnet, func$prodHTL),
          names.arg = c("New production", "Gross PP", "Net PP", "HTL"),
          xlab = TeX("Production (gC/m$^2$/yr)"),
          horiz=TRUE, las=1,
          border=NA)
}

simulate = function(p=parameters()) {
  #out = ode(c(0.1*p$N0, p$DOC0, p$POM0, p$B0), seq(0, p$tEnd, length.out = 100), derivative, p)
  out = cvode(seq(0, p$tEnd, length.out = p$tEnd),
              c(0.1*p$N0, p$DOC0, p$POM0, p$B0),
              function(t,y) derivative(t,y,p)[[1]],
              reltolerance = 1e-4,
              abstolerance = 1e-10+1e-4*c(0.1*p$N0, p$DOC0, p$POM0, p$B0))
  
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
    t = out[,1],
    y = out[,2:(p$n+4)],
    
    N = mean(out[ix,2]),
    DOC = mean(out[ix,3]),
    POM = mean(out[ix,4]),
    B = colMeans(out[ix,ixB]),
    
    Bmin = Bmin,
    Bmax = Bmax)
  
  result = c(result, list(rates = calcRates(max(result$t), result$N, result$DOC, result$POM, result$B,p)))
  return(result)
}

baserun = function(p = parameters()) {
  tic()
  sim = simulate(p)
  toc()
  
  plotSpectrum(sim)
  
  return(sim)
}
