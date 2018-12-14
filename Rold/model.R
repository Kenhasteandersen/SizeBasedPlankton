library(deSolve)
#
# Base parameters:
#
parameters <- function() {
  p = list()
  
  p$n = 20; # No of groups
  p$m = 10^seq(-8,1,length.out = p$n)  # Mass bins in mugC
  
  p$rhoCN = 5.68 # C:N ratio
  p$epsilonL = 0.8 # Light uptake efficiency
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
  
  # Camila uses 3.75e-5 (mugC)^-2/3 liter/(mugN day)
  p$AN = 0.005#0.000162 # Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  p$AL = 0.0019 # if using Al propto m^(2/3) for non-diatoms
  p$AL = 0.004864 # if using shading formula for non-diatoms
  p$cL = 0.08733 # if using shading formula for non-diatoms
  p$AF = 0.018  #  Fits to TK data for protists
  
  p$ANm = p$AN*p$m^(1/3)
  #p$ALm = p$AL*p$m^(2/3)*(1-nu)
  p$ALm = (1-nu) * p$cL*p$m * p$AL*p$m^(2/3) / ( p$cL*p$m + p$AL*p$m^(2/3) )  # shading formula
  p$AFm = (1-nu) * p$AF*p$m
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
  p$Jmax = 1.5 * p$m * (1-nu) # mugC/day
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
  p$L = 30  # PAR
  p$amplitudeL = 0 # amplitude of seasonal light variation in fractions of L
  
  p$tEnd = 365 # Simulation length (days)
  
  return(p)
}
#
# Food size preference:
#
phi = function(z, beta, sigma)
  exp( -(log(z/beta))^2/(2*sigma^2) )
#
# Seasonal light:
#
light = function(t, L, amplitudeL) {
  L*0.5*(1+amplitudeL*(-cos(t*2*pi/365))) 
}
#
# Seasonal diffusion:
#
diffusion = function(t, diff) {
  diff*(0.5*(1 + cos(2*t*pi/365)))^6
}
#
# Calc all rates in the cell model:
#
calcRates = function(t,N,DOC,POM,B,p) {
  with(p, {
    #
    # Potential uptakes:
    #
    JN =   ANm*N  # Diffusive nutrient uptake
    JDOC = ANm*DOC # Diffusive DOC uptake
    JL =   epsilonL*ALm*L*0.5*(1+amplitudeL*(-cos(t*2*pi/365)))  # Photoharvesting
    
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
    JLreal = JL-pmin(JCtot - (JF-JFreal)-f*Jmax, JL)
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
    
    diff = d
    if (amplitudeL>0)
      diff = diffusion(t, d)
    
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
      POMgeneration = POMgeneration,
      totKilled = sum(JFreal/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
  })
}
#
# Calc derivatives for the ode solver:
#
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
#
# Simulate the model:
#
simulate = function(p = parameters()) {
  # Simulate:
  nSave = 1000
  
  tEnd = p$tEnd
  if (p$amplitudeL>0)
    tEnd = 2*365
  
  out = ode(c(p$N0, p$DOC0, p$POM0, p$B0), 
            seq(0, tEnd, length.out = nSave), 
            derivative, p)
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
  
  result = c(result, 
             list(
               rates = calcRates(max(result$t), result$N, result$DOC, result$POM, result$B,p)))
  return(result)
}
#
# Calc ecosystem functions:
#
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
#==========Plots================

# Time-size plot of biomasses:
#
plotSizeTime = function(res) {
  t = res$t
  tlim = range(t)
  if (p$amplitudeL != 0) {
    t = res$t-max(res$t)+365
    tlim = c(0,365)
  }
  
  image(t, log10(res$p$m), log10(res$y[,4:13]), 
        xlim=tlim, zlim=c(-2,3),
        xlab="Time (days)",
        ylab="log10(cell mass)",
        col=colorRampPalette(c("blue","red"))(100))
}

plotResult = function(res, t=0) {
  p = res$p
  m = p$m
  
  #par(cex.axis=cex,
  #    cex.lab=cex,
  #    mar=c(4, 5, 6, 2) + 0.1)
  defaultplot(mfcol=c(2,1))
  
  alpha = 0.25
  colOsmo = rgb(0.5,0,0.5,alpha=alpha)
  colPhoto = rgb(0,1,0,alpha=alpha)
  colN = rgb(0,0,1,alpha=alpha)
  colMixo = rgb(1,0.5,0.5,alpha=alpha)
  colHetero = rgb(1,0,0,alpha=alpha)
  
  ylim = c(5,1000)
  fac = sqrt(m[2]/m[1])
  
  ixt = which(floor(res$t)==t+365)[1]
  if (p$amplitudeL==0) {
    L = p$L
    B = res$B
    N = res$N
    DOC = res$DOC
    r = res$rates
  } else {
    L = light(t, p$L, p$amplitudeL)
    B = res$y[ixt, 4:(p$n+3)]
    N = res$y[ixt, 1]
    DOC = res$y[ixt,2]
    r = calcRates(t, N, DOC, 0, B, p)
  }
  
  loglogpanel(ylim=ylim, xlim=m,
       xlab="Carbon mass ($\\mu$gC)",
       ylab="Biomass ($\\mu$gC/l)")
  if (p$amplitudeL==0)
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(res$Bmin, res$Bmax[seq(p$n,1,by = -1)]), 
            col=rgb(0.5,0.5,0.5,alpha=alpha), border=NA)
  lines(m, B, lwd=8)
  points(m,B)
  
  # Determine limiting process:
  ixOsmotroph = (r$JDOCreal > 5*r$JLreal)
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
  mtext(TeX("Diameter ($\\mu$m)"), side=3, line=3, at=1e-3, adj=1,cex=cex)
  
  legend(x="topright", bty="n", cex=cex,
         legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
         fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
         border=c("black","black","black","black","black","transparent"),
         lwd = c(0,0,0,0,0,3))
  
  text(x=m[1], y=10, labels=TeX(sprintf("DIN: %2.2f $\\mu$gN/l", N)) , cex=cex, pos=4, col=grey(0.5))
  text(x=m[1], y=7, labels=TeX(sprintf("DOC: %2.2f $\\mu$gC/l", DOC)), cex=cex, pos=4, col=grey(0.5))
  
  func = calcFunctions(res$p, res$rates, res$N, res$B)
  text(x=1e-2, 14, 
       labels=TeX(sprintf("Picoplankton: %2.2f $gC/m$^2$", func$Bpico)),
       cex=cex, pos=4, col=grey(0.5))
  text(x=1e-2, 10, 
       labels=TeX(sprintf("Nanoplankton: %2.2f $gC/m$^2$", func$Bnano)),
       cex=cex, pos=4, col=grey(0.5))
  text(x=1e-2, 7, 
       labels=TeX(sprintf("Microplankton: %2.2f $gC/m$^2$", func$Bmicro)),
       cex=cex, pos=4, col=grey(0.5))
  #text(x=m[1], y=7 , labels=sprintf("POM: %2.2f mugC/l",res$POM ), cex=cex, pos=4, col=grey(0.5))
  
  box()
  #
  # Rates
  #
  mm = 10^seq(-8,2, length.out = 100)  
 
  #par(cex.axis=cex,
  #    cex.lab=cex,
  #    mar=c(4, 5, 6, 2) + 0.1)
  ylim = c(0,1.5)

  semilogxpanel(xlim=m, ylim=ylim,
       xlab="Carbon mass ($\\mu$gC)",
       ylab="Rates (1/day)")
  lines(p$m, p$Jmax/p$m*r$f, lwd=8)
  lines(p$m, p$Jmax/p$m, lty=3)
  
  #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
  lines(p$m, p$ALm*L/p$m, lty=3, lwd=1, col="green")
  lines(p$m, r$JLreal/p$m, lwd=4, col="green")
  
  lines(mm, p$AN*mm^(1/3)*p$rhoCN*N/mm, lwd=1, lty=3, col="blue")
  lines(p$m, r$JNreal/p$m*p$rhoCN, lwd=4, col="blue")
  
  lines(mm, p$AN*mm^(1/3)*DOC/mm, lwd=1, lty=3, col="brown")
  lines(p$m, r$JDOCreal/p$m, lwd=4, col="brown")
  
  lines(p$m, r$JFreal/p$m,lwd=4,col="red")
  
  legend(x="topright", cex=cex,
         legend=c("Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
         col=c("green","blue","brown","red","black"),
         lwd=4,
         bty="n")
  
  
}