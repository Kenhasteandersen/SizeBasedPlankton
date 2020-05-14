source("model.R")
source("basetools.R")
library("sundialr")
library(tictoc)

parametersChemostat = function(p=parameters()) {
  #
  # Biogeochemical model:
  #
  p$d = 0.05  # diffusion rate, m/day
  p$M = 20   # Thickness of the mixed layer, m
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
  
  rates = calcRates(L,N,DOC,B,p)
  #
  # System:
  #
  if (p$latitude>0)
    diff = SeasonalExchange(p$latitude, t)
  else
    diff = p$d
  
  dBdt = diff*(0-B) + (rates$Jtot/p$m  - 
                         (  rates$mort+ 
                              rates$mortpred + 
                              rates$mort2 + 
                              p$mortHTL*p$mortHTLm))*B
  dBdt[(B<1e-3) & (dBdt<0)] = 0 # Impose a minimum concentration even if it means loss of mass balance
  
  #mortloss = sum(B*(rates$mort2 + p$mortHTL*(p$m>=p$mHTL)))
  mortloss = sum(B*((1-p$remin2)*rates$mort2 + p$mortHTL*p$mortHTLm))
  dNdt   =  diff*(p$N0-N) -
    sum(rates$JN*B/p$m) +
    sum(rates$JNloss*B/p$m) +
    p$remin2*sum(rates$mort2*B)/p$rhoCN +
    p$remin*mortloss/p$rhoCN
  dDOCdt =  diff*(0-DOC) -
    sum(rates$JDOC*B/p$m) +
    sum(rates$JCloss*B/p$m) +
    p$remin2*sum(rates$mort2*B) +
    p$remin*mortloss
  
  if (TRUE==FALSE) {
    # Check of nutrient conservation; should be close to zero
    Nin = diff*(p$N0-N) + diff*sum(0-B)/p$rhoCN
    Nout = (1-p$remin) * mortloss / p$rhoCN 
    NcheckSystem = Nin - Nout - sum(dBdt)/p$rhoCN - dNdt
    print(NcheckSystem)
    
    # Check carbon balance; should be close to zero:
    Cin = diff*(0-DOC) + diff*sum(0-B) + sum(rates$JLreal/p$m*B/p$epsilonL)
    Cout = (1-p$remin) * mortloss + sum(p$Jresp/p$m*B)
    CcheckSystem = Cin - Cout - sum(dBdt) - dDOCdt
    print(CcheckSystem)
  }
  #print(sum(rates$JF/p$epsilonF*B/p$m) - sum(rates$mortpred*B))
  
  # Check. Expensive to evaluate, so commented out  
  #  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
  #    browser()
  
  return(c(dNdt, dDOCdt, dBdt))
}

derivativeC = function(t,y,p) {
  derivC = .C("derivativeChemostat", 
              L=as.double(p$L), T=as.double(p$T), as.double(p$d), 
              as.double(p$N0), y=y, dudt=dudt)
  
  return(derivC$dudt)
}

compareCandRmodel = function(p,N=p$N0,DOC=p$DOC0,B=p$B0) {
  #
  # Remember to do a "simulate" first to reload library and set parameters
  #
  y = c(N,DOC,B)
  dudtR = derivative(0,y,p)
  
  # Load library
  dyn.load("../Cpp/model.so")
  # Set parameters
  dummy = .C("setParameters", as.integer(p$n), 
             p$m, p$rhoCN, p$epsilonL, p$epsilonF,
             p$ANm, p$ALm, p$AFm, p$Jmax, p$JFmaxm,
             p$Jresp, p$Jloss_passive_m,
             p$theta,
             p$mort, p$mort2, p$mortHTL*p$mortHTLm, p$remin,
             p$remin2, p$cLeakage);
  dudtC = derivativeC(0,y,p)
  
  plot(p$m, dudtR[3:(p$n+2)], type="l", log="x", col="red")
  lines(p$m, dudtC[3:(p$n+2)], col="blue", lty=dashed)
  print(dudtR[1:2])
  print(dudtC[1:2])
  
  return( calcRates(p$L,N,DOC,B,p) )
}

simulateChemostat = function(p=parametersChemostat(), useC=FALSE) {
  
  # Get the version of sundialr:
  pkg = installed.packages(fields = "Built")

  if (useC) {
    # Load library
    dyn.load("../Cpp/model.so")
    # Set parameters
    dummy = .C("setParameters", as.integer(p$n), 
               p$m, p$rhoCN, p$epsilonL, p$epsilonF,
               p$ANm, p$ALm, p$AFm, p$Jmax, p$JFmaxm,
               p$Jresp, p$Jloss_passive_m,
               p$theta,
               p$mort, p$mort2, p$mortHTL*p$mortHTLm, p$remin,
               p$remin2, p$cLeakage);
    
    dudt = assign("dudt", rep(0,p$n+2), envir = .GlobalEnv) # Need a static global for speed
    if (pkg[pkg[,1]=="sundialr"][3]=="0.1.2")
      out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                  IC = c(0.1*p$N0, p$DOC0, p$B0),
                  input_function = function(t,y, dummy) derivativeC(t,y,p),
                  reltolerance = 1e-6,
                  abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
    else
      out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                  IC = c(0.1*p$N0, p$DOC0, p$B0),
                  input_function = function(t,y, dummy) derivativeC(t,y,p),
                  reltolerance = 1e-6,
                  Parameters = 0,
                  abstolerance = 1e-10+1e-6*c(0.1*p$N0, p$DOC0, p$B0))
  } 
  else
    if (pkg[pkg[,1]=="sundialr"][3]=="0.1.2")
      out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                  IC = c(0.1*p$N0, p$DOC0, p$B0),
                  input_function = function(t,y, dummy) derivative(t,y,p),
                  reltolerance = 1e-5,
                  abstolerance = 1e-10+1e-5*c(0.1*p$N0, 1, p$B0))
    else
      out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
                  IC = c(0.1*p$N0, p$DOC0, p$B0),
                  input_function = function(t,y, dummy) derivative(t,y,p),
                  reltolerance = 1e-5,
                  Parameters = 0, 
                  abstolerance = 1e-10+1e-5*c(0.1*p$N0, 1, p$B0))
    
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
    
    result = c(result, list(rates = calcRates(SeasonalLight(p, max(result$t)), result$N, result$DOC, result$B,p)))
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
    prodCnet = conversion * sum( pmax(0, r$JLreal-r$JR)*B/m )
    #prodCnet = conversion * sum( r$JLreal*(1 - r$JR/(r$JCtot+r$JR))*B/m )
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
    resp = conversion * sum( r$JR*B/m )
    #
    # Bacterial production:
    #
    prodBact = conversion * sum( pmin(pmax(0,r$JDOC-r$JR), r$JNtot)*B/m )
    #prodBact = conversion * sum( r$JDOC*(1-r$JR/(r$JCtot+r$JR))*B/m )
    #
    # Efficiencies:
    #
    effHTL = prodHTL/prodNew # CHECK: CORRECT uNitS?
    effBact = prodBact / prodCnet
    #
    # Losses
    #
    lossPassive = conversion * sum( r$Jloss_passive/m*B ) 
    lossPhotouptake = conversion * sum( r$JCloss_photouptake/m*B )
    lossFeeding = conversion * sum( r$JCloss_feeding/m*B )
    lossFeedingHTL = sum( (1-epsilonF)*prodHTL )
    lossTotalC = lossPassive+lossPhotouptake+lossFeeding+lossFeedingHTL
    lossTotalN = (lossPassive+lossFeeding+lossFeedingHTL)/rhoCN
    #
    # Biomasses:
    #
    conversion = M/10*100*1e-6  # Convert to gC/m2
    d = calcESD(m)
    Bpico = conversion * sum( B[d < 2] )
    Bnano = conversion * sum( B[d>=2 & d <20])
    Bmicro = conversion * sum( B[d>=20])
    #
    # Chl-a (gC/m2): Use rough conversion from Edwards et al (2015) that Chl-a propto alpha
    #
    Chl_per_l = sum( r$JLreal/(epsilonL*L)/m*B )
    Chl_per_m2 = conversion * Chl_per_l
    
    return(list(
      prodNew = prodNew,
      prodCgross = prodCgross,
      prodCnet = prodCnet,
      prodHTL = prodHTL,
      prodBact = prodBact,
      resp = resp,
      effHTL = effHTL,
      effBact = effBact,
      lossPassive = lossPassive,
      lossPhotouptake = lossPhotouptake,
      lossFeeding = lossFeeding,
      lossFeedingHTL = lossFeedingHTL,
      lossTotalC = lossTotalC,
      lossTotalN = lossTotalN,
      Bpico=Bpico, Bnano=Bnano, Bmicro=Bmicro,
      Chl_per_m2 = Chl_per_m2,
      Chl_per_l = Chl_per_l
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
  #
  # Get functions
  # 
  func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
  # Get the func value from the previous call:
  oldfunc = attr(plotFunctionsChemostat, "oldfunc")
  if (is.null(oldfunc))
    oldfunc = func
  attr(plotFunctionsChemostat, "oldfunc") <<- func
  
  par(mfcol=c(2,1), mar=c(5,12,1,2))
  #
  # Fluxes:
  #
  flux = c(func$prodNew, func$prodCgross, func$prodCnet, func$prodHTL)
  oldflux = c(oldfunc$prodNew, oldfunc$prodCgross, oldfunc$prodCnet, oldfunc$prodHTL)
  
  heights = matrix(c(flux, oldflux), nrow=2, byrow = TRUE)
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
  #
  # Efficiencies:
  #
  PP = func$prodCnet
  eff = c(func$lossPassive/PP, func$lossPhotouptake/PP,
          func$lossFeeding/PP, func$lossTotalC/PP)
  PP = oldfunc$prodCnet
  oldeff = c(oldfunc$lossPassive/PP, oldfunc$lossPhotouptake/PP,
             oldfunc$lossFeeding/PP, oldfunc$lossTotalC/PP)
  heights = matrix(c(eff, oldeff), nrow=2, byrow=TRUE)
  barplot(height=heights,
          names.arg = c("ePassiveloss", "ePhotoloss", "eFeedingloss", "eTotalloss"),
          xlab = TeX("Fraction of net PP"),
          beside=TRUE, col=c("black","grey"),
          horiz=TRUE, las=1,
          border=NA)
}

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


plotSpectrum <- function(sim, t=max(sim$t), bPlot=TRUE) {
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
    r = calcRates(SeasonalLight(p,t), N, DOC, B, sim$p)
  }
  
  if (bPlot)
    defaultplot(mar=c(2.1,2.3,2.1,0))
  loglogpanel(xlim=p$m, ylim=ylim, 
              xlab="Carbon mass ($\\mu$gC)",
              ylab="Biomass ($\\mu$gC/l)")
  
  lines(m, B, lwd=8)
  if (p$n<15)
    points(m,b)
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
  text(x=10, 4.5, 
       labels=TeX(sprintf("Chl-a: %2.2f $mgC/m$^2$", 1000*func$Chl_per_m2)),
       cex=cex, pos=2, col=grey(0.5))
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

plotRates = function(sim, p=sim$p, 
                     B=sim$B, N=sim$N, DOC=sim$DOC,
                     t=max(sim$t), bPlot=TRUE) {
  mm = 10^seq(-8,2, length.out = 100)  
  
  L = p$L
  if (p$latitude!=0) {
    ixt = which(floor(sim$t)==t+365)[1]
    L = SeasonalLight(p,t)
    B = sim$y[ixt, 3:(p$n+2)]
    N = sim$y[ixt, 1]
    DOC = sim$y[ixt,2]
  }
  r = calcRates(L, N, DOC, B, p)
  
  if (bPlot)
    defaultplot()
  ylim = c(-1.4,1.6)
  semilogxpanel(xlim=p$m, ylim=ylim,
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Rates (1/day)")
  #
  # Gains
  #
  lines(p$m, r$Jtot/p$m, lwd=10, type="l", col="black")# log="x", xlim=range(p$m),
  lines(p$m, r$Jmax/p$m, lty=3)
  
  #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
  #JLreal = r$Jtot - r$JF+p$Jresp-r$JDOC
  #lines(p$m, p$ALm*p$L/p$m, lty=dotted, lwd=1, col="green")
  lines(p$m, r$JL/p$m , lty=dotted, lwd=1, col="green")
  lines(p$m, r$JLreal/p$m, lwd=4, col="green")
  
  #lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim$N/mm, lwd=1, lty=3, col="blue")
  #lines(p$m, p$Jmax * p$ANm*N / (p$Jmax/p$rhoCN + p$ANm*N)/p$m, lwd=1, lty=dotted, col="blue")
  #lines(p$m, r$JN/p$m*p$rhoCN, lwd=4, col="blue")
  #lines(p$m, p$JN/p$m, lwd=1, lty=dotted, col="blue")
  lines(p$m, r$JN/p$m, lwd=4, col="blue")
  
  #lines(mm, p$AN*mm^(1/3)*DOC/mm, lwd=1, lty=3, col="brown")
  lines(p$m, r$JDOC/p$m, lwd=4, col="brown")
  
  lines(p$m, r$JFreal/p$m,lwd=4,col="red")
  lines(p$m, p$epsilonF * r$JFmax/p$m ,col="red", lty=dotted)
  
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
  
  JNexude = r$JNloss 
  lines(p$m, -(r$mortpred + p$mortHTL*p$mortHTLm + r$mort2 + p$mort), lwd=10)
  lines(p$m, -r$mortpred, col="red", lwd=4)
  lines(p$m, -p$mortHTL*p$mortHTLm, col="magenta", lwd=4)
  lines(p$m, -r$mort2, col="orange", lwd=4)
  lines(p$m, -r$JR/p$m, col="grey", lwd=4)
  lines(p$m, -r$Jloss_passive/p$m, col="darkgreen", lwd=4)
  
  BSheldon =exp(mean(log(B)))
  delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
  mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / delta
  lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")
  
  legend(x="bottomright", cex=cex,
         legend=c("Losses:", "Predation", "Virulysis", 
                  "Higher trophic levels","Respiration","Passive"),
         col=c(NA,"red", "orange", "magenta","grey","darkgreen"),
         lwd=c(0,4,4,4,4,4), bty="n")
  
  lines(p$m, 0*p$m, col="white", lwd=4)
  lines(p$m, 0*p$m, lty=dashed, lwd=2)
  
  return(r)
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
    r = calcRates(L,N, DOC, B, p)
  }
  
  defaultplot()
  semilogxpanel(xlim=m, ylim=c(0,0.4),
                xlab="Carbon mass ($\\mu$gC)",
                ylab="Loss rates (1/day)")
  
  lines(m, (r$JNloss*p$rhoCN-r$JCloss_feeding)/m, col="blue", lwd=4)
  #lines(m, (r$JNloss-r$JNloss_piss)/m, col="blue", lwd=3, lty=dashed)
  lines(m, r$JCloss_feeding/m, col="red", lwd=4)
  lines(m, r$JCloss_photouptake/m, col="green", lwd=4)
  lines(m, r$Jloss_passive/m, col="darkgreen", lwd=4)
  
  legend(x="topright", cex=cex,
         legend=c("Leaks:","N exudation", "C exudation",
                  "N+C sloppy feeding","Passive exudation"),
         col=c("white","blue","green","red","darkgreen"),
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
    r = calcRates(L, N, DOC, B, p)
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
  lines(m, -r$JR/m, col="magenta")
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
  sim = simulateChemostat(p, useC)
  toc()
  
  defaultplot(c(2,1))
  plotSpectrum(sim, bPlot=FALSE)
  plotRates(sim, bPlot=FALSE)
  
  return(sim)
}
