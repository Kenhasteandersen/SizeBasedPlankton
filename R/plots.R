require(latex2exp)
source("basetools.R")
source("model.R")

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
  
  ylim = c(1,500)
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
  #
  # Add gray-scale variation
  #
  if (p$amplitudeL==0)
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim$Bmin, sim$Bmax[seq(p$n,1,by = -1)]), 
            col=rgb(0.5,0.5,0.5,alpha=alpha), border=NA)
  
  # Determine limiting process:
  ixOsmotroph = (r$JDOCreal > r$JLreal)
  ixLlimited = ((r$JCtot < r$JNtot*p$rhoCN) & !ixOsmotroph)
  ixHetero = (r$JNloss_piss>0)
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
         lwd = c(0,0,0,0,0,3),
         col=c(NA,NA,NA,NA,NA,1))
  #
  # Summary state variables: 
  #
  text(x=m[1], y=2.25, labels=TeX(sprintf("DIN: %2.2f $\\mu$mol N/l", N/14)) , cex=cex, pos=4, col=grey(0.5))
  text(x=m[1], y=1.5, labels=TeX(sprintf("DOC: %2.2f $mmol C/l", 1000*DOC/12)), cex=cex, pos=4, col=grey(0.5))
  
  func = calcFunctions(sim$p, sim$rates, sim$N, sim$B)
  text(x=10, 3.3, 
       labels=TeX(sprintf("Picoplankton: %2.2f $mgC/m$^2$", 1000*func$Bpico)),
       cex=cex, pos=2, col=grey(0.5))
  text(x=10, 2.25, 
       labels=TeX(sprintf("Nanoplankton: %2.2f $mgC/m$^2$", 1000*func$Bnano)),
       cex=cex, pos=2, col=grey(0.5))
  text(x=10, 1.5, 
       labels=TeX(sprintf("Microplankton: %2.2f $mgC/m$^2$", 1000*func$Bmicro)),
       cex=cex, pos=2, col=grey(0.5))

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
  ylim = c(-1.4,1.6)
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
         legend=c("Gains:","Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
         col=c("white","green","blue","brown","red","black"),
         lwd=c(0,4,4,4,4,4),
         bty="n")
  #
  # Mortality losses
  #
  polygon(c(1e-9,10,10,1e-9), c(-1.5,-1.5,0,0), 
          col=rgb(1,0,0,alpha=0.25), border=NA)
  
  mortHTL = p$mortHTL*(p$m>p$mHTL)
  
  lines(p$m, -(sim$rates$mortpred + mortHTL + sim$rates$mort2 + p$mort), lwd=10)
  lines(p$m, -sim$rates$mortpred, col="red", lwd=4)
  lines(p$m[p$m>=p$mHTL], -p$mortHTL*sign(p$m[p$m>p$mHTL]), col="magenta", lwd=4)
  lines(p$m, -sim$rates$mort2, col="orange", lwd=4)
  
  BSheldon =exp(mean(log(sim$B)))
  delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
  mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / delta
  lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")
  
  legend(x="bottomright", cex=cex,
         legend=c("Losses:", "Predation", "Virulysis", "Higher trophic levels"),
         col=c(NA,"red", "orange", "magenta"),
         lwd=c(0,4,4,4), bty="n")
  
  lines(p$m, 0*p$m, col="white", lwd=4)
  lines(p$m, 0*p$m, lty=dashed, lwd=2)
}

plotLeaks = function(sim, t=max(sim$t)) {
  p = sim$p
  
  m = p$m
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
  semilogxpanel(xlim=m, ylim=c(0,0.4),
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Loss rates (1/day)")
  
  lines(m, r$JNloss_piss*p$rhoCN/m, col="blue", lwd=4)
  #lines(m, (r$JNloss-r$JNloss_piss)/m, col="blue", lwd=3, lty=dashed)
  lines(m, r$JCloss_feeding/m, col="red", lwd=4)
  lines(m, (r$JCloss-r$JCloss_feeding)/m, col="green", lwd=4)
  
  legend(x="topright", cex=cex,
         legend=c("Leaks:","N exudation", "C exudation", "N+C sloppy feeding"),
         col=c("white","blue","green","red"),
         lwd=4, bty="n")
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