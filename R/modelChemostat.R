source("model.R")
library("sundialr")

parametersChemostat = function(p=parameters()) {
  #
  # Biogeochemical model:
  #
  p$d = 0.05  # diffusion rate, m/day
  p$M = 30   # Thickness of the mixed layer, m
  p$T = 10   # Temperature
  p$N0 = 150 # Deep nutrient levels
  
  #
  # Light:
  #
  p$L = 60  # PAR, mu E/m2/s
  p$latitude = 0 # amplitude of seasonal light variation in fractions of L
  
  p$tEnd = 365 # Simulation length (days)
  
  return(p)
}

#
# Seasonal variation in exchange rate as a function of latitude (degrees)
# and time (days)
#
SeasonalExchange = function(latitude, t) {
  t = t %% 365
  
  dmax = 0.05*(1+tanh(0.05*(latitude-40)))
  dsummer = 0.01
  tspring = 180 * latitude/120
  tautumn = 200 + 180 *(90-latitude)/90
  widthautumn = 1
  
  summer = 1-0.5*(1+tanh(8*((t-tspring)/365*2*pi)))
  winter = 1-0.5*(1 - tanh(widthautumn*((t-tautumn)/365*2*pi)))
  spring = 1-0.5*(1 - tanh(widthautumn*(((t+365)-tautumn)/365*2*pi)))
  summer = pmin(summer, spring)
  
  d = dsummer + dmax*(winter + summer)
}
#
# Seasonal variation in light. Roughly taken from Evans and Parslows. 
# M is the depth of the mixed layer.
#
SeasonalLight = function(p,t) {
  p$L*exp(-0.025*p$M)*(1 - 0.8*sin(pi*p$latitude/180)*cos(2*pi*t/365))
}

derivative = function(t,y,p) {
  N = y[1]
  DOC = y[2]
  B = y[3:(2+p$n)]
  
  if (p$latitude > 0)
    L = SeasonalLight(p,t)
  else
    L = p$L
  
  rates = calcRates(t,L,N,DOC,B,p)
  #
  # System:
  #
  if (p$latitude>0)
    diff = SeasonalExchange(p$latitude, t)
  else
    diff = p$d
  
  dBdt = diff*(0-B) + (rates$Jtot/p$m  - 
                         (rates$mort+ rates$mortpred + rates$mort2 + p$mortHTL*(p$m>=p$mHTL)))*B
  dBdt[(B<1e-3) & (dBdt<0)] = 0 # Impose a minimum concentration even if it means loss of mass balance
  mortloss = sum(B*(rates$mort2 + p$mortHTL*(p$m>=p$mHTL)))
  dNdt   =  diff*(p$N0-N)  - sum(rates$JN*B/p$m)   + sum(rates$JNloss*B/p$m) +
    p$remin*mortloss/p$rhoCN
  dDOCdt =  diff*(0-DOC) - sum(rates$JDOC*B/p$m) + sum(rates$JCloss*B/p$m) + 
    p$remin*mortloss
  
  # Check. Expensive to evaluate, so commented out  
  #  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
  #    browser()
  
  return(c(dNdt, dDOCdt, dBdt))
}

dudt = rep(0,12) # Need a static global for speed
derivativeC = function(t,y,p) {
  derivC = .C("derivativeChemostat", 
              L=as.double(p$L), T=as.double(p$T), as.double(p$d), 
              as.double(p$N0), y=y, dudt=dudt)
  
  return(derivC$dudt)
}

simulate = function(p=parametersChemostat(), useC=FALSE) {
  if (useC) {
    # Load library
    dyn.load("../Cpp/model.so")
    # Set parameters
    dummy = .C("setParameters", as.integer(p$n), 
               p$m, p$rhoCN, p$epsilonL, p$epsilonF,
               p$ANm, p$ALm, p$AFm, p$Jmax, p$Jresp, p$theta,
               p$mort, p$mort2, 0*p$m + p$mortHTL*(p$m>p$mHTL), p$remin);
    
    out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                IC = c(0.1*p$N0, p$DOC0, p$B0),
                input_function = function(t,y) derivativeC(t,y,p),
                reltolerance = 1e-6,
                abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  } else
    out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                IC = c(0.1*p$N0, p$DOC0, p$B0),
                input_function = function(t,y) derivative(t,y,p),
                reltolerance = 1e-6,
                abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  
  nSave = dim(out)[1]
  # Assemble results:
  ix = seq(floor(nSave/2),nSave)
  ixB = 4:(p$n+3)
  
  Bmin = 0*p$m
  Bmax = 0*p$m
  for (i in 1:p$n) {
    Bmin[i] = max(1e-20, min(out[ix,ixB[i]]))
    Bmax[i] = max(1e-20, out[ix,ixB[i]])
  }
  
  result = list(
    p = p,
    t = out[,1],
    y = out[,2:(p$n+3)],
    
    N = mean(out[ix,2]),
    DOC = mean(out[ix,3]),
    B = colMeans(out[ix,ixB]),
    
    Bmin = Bmin,
    Bmax = Bmax)
  
  result = c(result, list(rates = calcRates(max(result$t), p$L, result$N, result$DOC, result$B,p)))
  return(result)
}


calcFunctionsChemostat = function(param,r,N,B) {
  with(param, {
    conversion  = 365*M*1e-6*1000 # Convert to gC/yr/m2
    #
    # New production calculated as the flux of nitrogen into the system. Units of carbon:
    #
    prodNew = conversion * d*(N0-N) * rhoCN  
    #
    # Primary production (carbon fixed)
    #
    prodCgross = conversion * sum(r$JLreal*B/m)/epsilonL
    prodCnet = conversion * sum( (r$JLreal-Jresp)*B/m )
    #
    # Loss to HTL:
    #
    prodHTL = conversion * sum( mortHTL*(m>=mHTL)*B )
    #
    # Loss to depth:
    #
    prodSeq = 0#conversion * r$POMgeneration * (1-epsilonPOM)
    #
    # Respiration
    #
    resp = conversion * sum( r$Jresp*B/m )
    #
    # Efficiencies:
    #
    effHTL = prodHTL/prodNew # CHECK: CORRECT uNitS?
    #
    # Biomasses:
    #
    conversion = M/10*100*1e-6  # Convert to gC/m2
    d = calcESD(m)
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

plotTimeline = function(sim, time=max(sim$t)) {
  p = sim$p
  t = sim$t
  
  par(cex.axis=cex,
      cex.lab=cex,
      mar=c(4, 5, 6, 2) + 0.1)
  
  y = sim$y
  y[y <= 0] = 1e-30
  
  ylim = c(max(1e-5, min(sim$y)), max(sim$y))
  if (p$latitude==0) {
    xlim = range(t)  
  } else {
    xlim = c(0,365)
    t = 365+t-max(t)
  }
  
  plot(t, y[,1], log="y", type="l", col="blue", 
       ylim=ylim, xlim=xlim, lwd=2,
       xlab="Time (day)", ylab=TeX("Biomass ($\\mu$gC/l)"))
  lines(t, y[,2], col="magenta", lwd=2)
  lines(t, y[,3], col="orange", lwd=2)
  for (i in 1:p$n)
    lines(t, y[,i+2], lwd=i/p$n*3, col="black")
  
  if (p$latitude>0) {
    lines(time*c(1,1), ylim, lty=dotted)
    lines(t, SeasonalLight(p,t),
          col="orange", lwd=2)
  }
}

plotSeasonal = function(p,time) {
  defaultplot()
  defaultpanel(xlim=c(0,365),
               ylim=c(0,1),
               xlab="Time (days)",
               ylab="d/max(d) / L/max(L)")
  
  t = seq(0,365,length.out = 100)
  lines(t, SeasonalExchange(p$latitude, t)/max(SeasonalExchange(p$latitude, t)),col="red", lwd=3)
  lines(t, SeasonalLight(p, t)/p$L, col="green", lwd=3)
  vline(x=time)
}

plotSeasonalTimeline = function(sim) {
  require(lattice)
  require(latticeExtra)
  
  p = sim$p
  ix = sim$t>max(sim$t-365)
  t = sim$t[ix]
  
  B = log10(sim$y[ix,3:(p$n+2)])
  B[B<(-1)] = -1
  B[B>3] = 3
  
  levelplot(
    B, 
    column.values=(p$m), 
    row.values = t, 
    aspect="fill",
    scales = list(y = list(log = 10)),
    yscale.components = yscale.components.log10ticks,
    xlab = "Time (days)",
    ylab = TeX("Carbon mass ($\\mu$gC)"),
    col.regions = terrain.colors(100),
    ylim=c(p$m[1]/3, max(p$m)*3)
  )
}

#
# Plot functions:
#
plotFunctionsChemostat <- function(sim) {
  # Get the func value from the previous call:
  oldfunc = attr(plotFunctionsChemostat, "oldfunc")
  if (is.null(oldfunc))
    oldfunc = c(0,0,0,0)
  
  func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
  fnc = c(func$prodNew, func$prodCgross, func$prodCnet, func$prodHTL)
  attr(plotFunctionsChemostat, "oldfunc") <<- fnc
  heights = matrix(c(fnc, oldfunc), nrow=2, byrow = TRUE)
  
  par(mar=c(5,12,4,2))
  barplot(height=heights,
          names.arg = c("New production", "Gross PP", "Net PP", "HTL"),
          xlab = TeX("Production (gC/m$^2$/yr)"),
          beside=TRUE, col=c("black","grey"),
          horiz=TRUE, las=1,
          border=NA)
  legend("topright",
         c("This simulation","Previous simulation"),
         fill=c("black","grey"),
         bty="n")
}

source("basetools.R")
source("model.R")


#
# Returns the trophic strategy as one of: osmoheterotroph, light or nutrient-limited
# phototroph, mixotroph, or heterotroph.
#
calcStrategy = function(p,r) {
  strategy = rep('Unknown', p$n)
  strategy[r$JN*p$rhoCN>r$JL] = "Light limited"
  strategy[r$JL>=r$JN*p$rhoCN] = "Nutrient limited"
  strategy[r$JDOC > r$JL] = "Osmoheterotroph"
  strategy[(r$JNloss>1e-5) & (r$JN<r$JF/p$rhoCN)] = "Heterotroph"
  strategy[(r$JF/r$JL > 0.25) & !strategy=="Heterotroph"] = "Mixotroph"  
  return(strategy)
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
  if (p$latitude==0) {
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, B, p)
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
  if (p$latitude==0)
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim$Bmin, sim$Bmax[seq(p$n,1,by = -1)]), 
            col=rgb(0.5,0.5,0.5,alpha=alpha), border=NA)
  
  # Determine limiting process:
  strategy = calcStrategy(p,r)
  for (i in 1:p$n) {
    if (strategy[i]=="Heterotroph")
      col = colHetero
    if (strategy[i]=="Mixotroph")
      col = colMixo
    if (strategy[i]=="Light limited")
      col = colPhoto
    if (strategy[i]=="Nutrient limited")
      col = colN
    if (strategy[i]=="Osmoheterotroph")
      col = colOsmo
    
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
  
  func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
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
  if (p$latitude==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, B, p)
  }
  
  defaultplot()
  ylim = c(-1.4,1.6)
  semilogxpanel(xlim=p$m, ylim=ylim,
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Rates (1/day)")
  #
  # Gains
  #
  lines(p$m, r$Jtot/p$m, lwd=10, type="l", col="black")# log="x", xlim=range(p$m),
  lines(p$m, p$Jmax/p$m, lty=3)
  
  #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
  #JLreal = r$Jtot - r$JF+p$Jresp-r$JDOC
  #lines(p$m, p$ALm*p$L/p$m, lty=dotted, lwd=1, col="green")
  lines(p$m, p$epsilonL*p$Jmax*p$ALm*L / (p$Jmax + p$ALm*L)/p$m , lty=dotted, lwd=1, col="green")
  lines(p$m, r$JLreal/p$m, lwd=4, col="green")
  
  #lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim$N/mm, lwd=1, lty=3, col="blue")
  lines(p$m, p$Jmax * p$ANm*N / (p$Jmax/p$rhoCN + p$ANm*N)/p$m, lwd=1, lty=3, col="blue")
  lines(p$m, r$JN/p$m*p$rhoCN, lwd=4, col="blue")
  
  lines(mm, p$AN*mm^(1/3)*sim$DOC/mm, lwd=1, lty=3, col="brown")
  lines(p$m, r$JDOC/p$m, lwd=4, col="brown")
  
  lines(p$m, r$JF/p$m,lwd=4,col="red")
  
  legend(x="topright", cex=cex,
         legend=c("Gains:","Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
         col=c("white","green","blue","brown","red","black"),
         lwd=c(0,4,4,4,4,4),
         bty="n")
  #
  # Losses
  #
  polygon(c(1e-9,10,10,1e-9), c(-1.5,-1.5,0,0), 
          col=rgb(1,0,0,alpha=0.25), border=NA)
  
  mortHTL = p$mortHTL*(p$m>p$mHTL)
  JNexude = r$JNloss 
  lines(p$m, -(sim$rates$mortpred + mortHTL + sim$rates$mort2 + p$mort), lwd=10)
  lines(p$m, -sim$rates$mortpred, col="red", lwd=4)
  lines(p$m[p$m>=p$mHTL], -p$mortHTL*sign(p$m[p$m>p$mHTL]), col="magenta", lwd=4)
  lines(p$m, -sim$rates$mort2, col="orange", lwd=4)
  lines(p$m, -p$Jresp/p$m, col="grey", lwd=4)
  
  BSheldon =exp(mean(log(sim$B)))
  delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
  mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / delta
  lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")
  
  legend(x="bottomright", cex=cex,
         legend=c("Losses:", "Predation", "Virulysis", "Higher trophic levels","Respiration"),
         col=c(NA,"red", "orange", "magenta","grey"),
         lwd=c(0,4,4,4,4), bty="n")
  
  lines(p$m, 0*p$m, col="white", lwd=4)
  lines(p$m, 0*p$m, lty=dashed, lwd=2)
}

plotLeaks = function(sim, t=max(sim$t)) {
  p = sim$p
  
  m = p$m
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$latitude==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, B, p)
  }
  
  defaultplot()
  semilogxpanel(xlim=m, ylim=c(0,0.4),
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Loss rates (1/day)")
  
  lines(m, (r$JNloss*p$rhoCN-r$JCloss_feeding)/m, col="blue", lwd=4)
  #lines(m, (r$JNloss-r$JNloss_piss)/m, col="blue", lwd=3, lty=dashed)
  lines(m, r$JCloss_feeding/m, col="red", lwd=4)
  lines(m, r$JCloss_photouptake/m, col="green", lwd=4)
  
  legend(x="topright", cex=cex,
         legend=c("Leaks:","N exudation", "C exudation", "N+C sloppy feeding"),
         col=c("white","blue","green","red"),
         lwd=4, bty="n")
}

plotComplexRates = function(sim, t=max(sim$t)) {
  p = sim$p
  
  ixt = which(floor(sim$t)==t+365)[1]
  if (p$latitude==0) {
    L = p$L
    B = sim$B
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
  } else {
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
    r = calcRates(t, N, DOC, B, p)
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


baserunChemostat = function(p = parametersChemostat(), useC=FALSE) {
  tic()
  sim = simulate(p, useC)
  toc()
  plotSpectrum(sim)
  return(sim)
}
