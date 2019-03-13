library("sundialr")
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
  p$ALm = p$AL*p$m^(2/3) * (1-exp(- p$cL*p$m^(1/3) ))  # shading formula
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
  p$alphaJ = 1.5
  p$Jmax = p$alphaJ * p$m * (1-nu) # mugC/day
  p$cR = 0.1
  p$Jresp = p$cR*p$Jmax
  #
  # Losses:
  #
  p$mort = 0*0.005*(p$Jmax/p$m) * p$m^(-1/4)
  p$mort2 = 0.0002*p$n
  p$mortHTL = 0.2
  p$mHTL = max(p$m)/p$beta # Bins affecte by HTL mortality
  
  p$remin = 0.0 # fraction of mortality losses reminerilized to N and DOC
  #
  # Biogeochemical model:
  #
  p$d = 0.05  # diffusion rate, m/day
  p$M = 30   # Thickness of the mixed layer, m
  p$N0 = 150 # Deep nutrient levels
  #
  # Initial conditions:
  #
  p$DOC0 = 0
  p$B0 = rep(10,p$n)
  #
  # Light:
  #
  p$L = 20  # PAR, mu E/m2/s
  p$amplitudeL = 0 # amplitude of seasonal light variation in fractions of L
  
  p$tEnd = 365 # Simulation length (days)
  
  return(p)
}

phi = function(z, beta, sigma)
  exp( -(log(z/beta))^2/(2*sigma^2) )

calcRates = function(t,N,DOC,B,p) {
  with(p, {
    B = pmax(0,B)
    #
    # Uptakes
    #
    JN =   Jmax/p$rhoCN * ANm*N / (Jmax/p$rhoCN + ANm*N) # Diffusive nutrient uptake
                                                        # in units of N/time
    JDOC = Jmax * ANm*DOC / (Jmax + ANm*DOC) # Diffusive DOC uptake, units of C/time
    
    LL = L*(1+amplitudeL*(-cos(t*2*pi/365)))
    JL =   epsilonL * Jmax * ALm*LL / (Jmax + ALm*LL)  # Photoharvesting
    
    F = theta %*% B
    JF = epsilonF * Jmax * AFm*F / (Jmax + AFm*F)        # Feeding
    
    JCtot = JL+JF-Jresp+JDOC # Total carbon intake
    JNtot = JN+JF/rhoCN # In units of N
    Jtot = pmin( JCtot, JNtot*rhoCN )  # Liebigs law; units of C
    
    JLreal = Jtot - JF+p$Jresp-JDOC
    # 
    # Losses:
    #
    Jloss_feeding = (1-epsilonF)/epsilonF*JF # Incomplete feeding (units of carbon per time)
    JNlossLiebig = pmax(0, JNtot*rhoCN-JCtot)/rhoCN  # N losses from Liebig
    JClossLiebig = pmax(0, JF-Jresp+JDOC-JNtot*rhoCN) # C losses from Liebig, not counting losses from photoharvesting
    #JClossLiebig = pmin(JClossLiebig, JDOC) # However, light surplus is not leaked but is downregulated

    JNloss = Jloss_feeding/rhoCN + JNlossLiebig
    JCloss = Jloss_feeding + JClossLiebig + (1-epsilonL)/epsilonL*JLreal
    
    #if (sum(c(JNloss,JCloss,B)<0))
    #  browser()
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JF/epsilonF*B/m/F)
    #
    # System:
    #
    
    #dSeason = function(t,dmax) 
    #  0.5*(1+tanh(8*(2-t/365*2*pi))) + 0.5*(1 - tanh(2*(5-t/365*2*pi)))
    #d = d + dSeason(t,1)
    
    dBdt = (Jtot/m  - (mort+ mortpred + mort2*B + mortHTL*(m>=mHTL)))*B
    mortloss = sum(B*(mort2*B + mortHTL*(m>=mHTL)))
    dNdt   =  d*(N0-N)  - sum(JN*B/m)   + sum(JNloss*B/m) + remin*mortloss/rhoCN
    dDOCdt =  d*(0-DOC) - sum(JDOC*B/m) + sum(JCloss*B/m) + remin*mortloss
    
    # Check of nutrient conservation; should be close to zero
    #Nin = d*(N0-N)
    #Nout = (1-remin) * mortloss / rhoCN
    #NcheckSystem = Nin - Nout - sum(dBdt)/rhoCN - dNdt
    #print(NcheckSystem)
    
    #print(sum(JF/epsilonF*B/m) - sum(mortpred*B))
    
    return(list( 
      dNdt=dNdt, dDOCdt=dDOCdt, dBdt=dBdt, 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JLreal = JLreal,
      JNlossLiebig=JNlossLiebig, JClossLiebig=JClossLiebig,
      Jloss_feeding=Jloss_feeding, JCloss=JCloss, JNloss=JNloss,
      Jtot=Jtot, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      totKilled = sum(JF/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jtot*B/m)))  
  })
}
#
# Version 1. Heuristic down-regulation; first of feeding,
#
calcRatesOld = function(t,N,DOC,B,p) {
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
    
    f = Jtot / (Jtot + Jmax) # feeding level
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
    JNloss_piss = -pmin(0, (f*Jmax - JFreal))/rhoCN  # Exudation of surplus nutrients by heterotrophs
    JNloss = (1-epsilonF)/epsilonF*JFreal/rhoCN + JNloss_piss
    JCloss_feeding = (1-epsilonF)/epsilonF*JFreal
    JCloss = JCloss_feeding + (1-epsilonL)/epsilonL*JLreal
    
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
    dBdt = (Jmax*f/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL)))*B
    dNdt   =  d*(N0-N)  - sum(JNreal*B/m)   + sum(JNloss*B/m)
    dDOCdt =  d*(0-DOC) - sum(JDOCreal*B/m) + sum(JCloss*B/m)

    # Check of nutrient conservation; should be close to zero
    #Nin = d*(N0-N)
    #Nout = (1-remin) * ( sum(mort2*B*B) + sum(mortHTL*(m>=mHTL)*B) ) / rhoCN
    #NcheckSystem = Nin - Nout - sum(dBdt)/rhoCN - dNdt
    #print(NcheckSystem)
    
    return(list( 
      dNdt=dNdt, dDOCdt=dDOCdt, dBdt=dBdt, 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JNreal=JNreal, JDOCreal=JDOCreal, JLreal=JLreal, JFreal=JFreal, 
      JNloss_piss=JNloss_piss, JNloss=JNloss, 
      JCloss_feeding=JCloss_feeding, JCloss=JCloss,
      Jtot=Jtot, f=f, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      totKilled = sum(JFreal/epsilonF*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
  })
}

derivative = function(t,y,p) {
  N = y[1]
  DOC = y[2]
  B = y[3:(2+p$n)]
  
  rates = calcRates(t,N,DOC,B,p)
  
  # Check. Expensive to evaluate, so commented out  
  #  if ( sum(c( is.nan(unlist(rates)), is.infinite(unlist(rates)), rates$N<0, rates$B<0))>0)
  #    browser()
  
  return(list(c(rates$dNdt, rates$dDOCdt, rates$dBdt)))
}

simulate = function(p=parameters()) {
  out = cvode(time_vector = seq(0, p$tEnd, length.out = p$tEnd),
              IC = c(0.1*p$N0, p$DOC0, p$B0),
              input_function = function(t,y) derivative(t,y,p)[[1]],
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
  
  result = c(result, list(rates = calcRates(max(result$t), result$N, result$DOC, result$B,p)))
  return(result)
}

baserun = function(p = parameters()) {
  tic()
  sim = simulate(p)
  toc()
  
  return(sim)
}
