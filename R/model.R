

# --------------------------------------------------
# Core logic for the model
# --------------------------------------------------

parameters <- function() {
  p = list()
  
  p$n = as.integer(10); # No of groups
  p$m = 10^seq(-8,1,length.out = p$n)  # Mass bins in mugC
  
  p$rhoCN = 5.68 # C:N mass ratio
  p$epsilonL = 0.9 # Light uptake efficiency
  p$epsilonF = 0.8 # Assimilation efficiency
  p$cLeakage = 0.00015 # passive leakage of C and N
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
  p$alphaJ = 1.5 # per day
  p$Jmax = p$alphaJ * p$m * (1-nu) # mugC/day
  p$cR = 0.1
  p$Jresp = p$cR*p$alphaJ*p$m
  #
  # Losses:
  #
  p$mort = 0*0.005*(p$Jmax/p$m) * p$m^(-1/4)
  p$mort2 = 0.0002*p$n
  p$mortHTL = 0.1
  p$mHTL = max(p$m)/p$beta # Bins affected by HTL mortality
  
  p$remin = 0.0 # fraction of mortality losses reminerilized to N and DOC
  p$remin2 = 0.3 # fraction of virulisus remineralized to N and DOC
  
  p$T = 10
  #
  # Initial conditions:
  #
  p$N0 = 150
  p$DOC0 = 0
  p$B0 = rep(1,p$n)
  
  return(p)
}
#
# Prey size preference function:
#
phi = function(z, beta, sigma) {
  exp( -(log(z/beta))^2/(2*sigma^2) )
}
#
# Q10 temperature function:
#
fTemp = function(Q10, T) {
  return(Q10^(T/10-1))
}
#
# Convert to ESD:
#
calcESD = function(m) {
  10000 * 1.5 * (m*1e-6)^(1/3)
}

calcRates = function(t,L,N,DOC,B,p) {
  with(p, {
    B = pmax(0,B)
    N = max(0,N)
    DOC = max(0, DOC)
    #
    # Temperature corrections:
    #
    ANmT = ANm*fTemp(1.5,p$T)
    JmaxT = Jmax*fTemp(2,p$T)
    JR = Jresp*fTemp(2,p$T)
    #
    # Uptakes
    #
    JN =   JmaxT/p$rhoCN * ANmT*N / (JmaxT/p$rhoCN + ANmT*N) # Diffusive nutrient uptake
                                                        # in units of N/time
    JDOC = JmaxT * ANmT*DOC / (JmaxT + ANmT*DOC) # Diffusive DOC uptake, units of C/time
    
    JL =   epsilonL * JmaxT * ALm*L / (JmaxT + ALm*L)  # Photoharvesting
    
    F = theta %*% B
    JF = epsilonF * JmaxT * AFm*F / (JmaxT + AFm*F)        # Feeding

    # Total nitrogen uptake:    
    JNtot = JN+JF/rhoCN # In units of N
    
    # Down-regulation of light uptake:
    JLreal = pmin(JL, pmax(0, JNtot*rhoCN - (JDOC-JR)))
                     
    JCtot = JLreal+JF+JDOC-JR # Total carbon uptake

    Jtot = pmin( JCtot, JNtot*rhoCN )  # Liebigs law; units of C
    # 
    # Losses:
    #
    JCloss_feeding = (1-epsilonF)/epsilonF*JF # Incomplete feeding (units of carbon per time)
    JCloss_photouptake = (1-epsilonL)/epsilonL*JLreal
    JNlossLiebig = pmax(0, JNtot*rhoCN-JCtot)/rhoCN  # N losses from Liebig
    JClossLiebig = pmax(0, JCtot-JNtot*rhoCN) # C losses from Liebig, not counting losses from photoharvesting
    #JClossLiebig = pmax(0, Jtot - JNtot*rhoCN) # C losses from Liebig, not counting losses from photoharvesting
    #JClossLiebig = pmin(JClossLiebig, JDOC) # However, light surplus is not leaked but is downregulated

    Jloss_passive = p$cLeakage * m^(2/3) # in units of C
    
    JNloss = JCloss_feeding/rhoCN + JNlossLiebig + Jloss_passive/rhoCN
    JCloss = JCloss_feeding + JCloss_photouptake + JClossLiebig + Jloss_passive
    
    #if (sum(c(JNloss,JCloss,B)<0))
    #  browser()
    #
    # Mortality:
    #
    mortpred =  t(theta) %*% (JF/epsilonF*B/m/F)

    return(list( 
      JN=JN, JDOC=JDOC, JL=JL, JF=JF,
      JLreal = JLreal, JR=JR,
      JNlossLiebig=JNlossLiebig, JClossLiebig=JClossLiebig,
      JCloss_photouptake=JCloss_photouptake,
      Jloss_passive = Jloss_passive,
      JCloss_feeding=JCloss_feeding, JCloss=JCloss, JNloss=JNloss,
      Jtot=Jtot, F=F, JCtot = JCtot, JNtot=JNtot,
      mortpred=mortpred, mort=mort,
      mort2=mort2*B,
      totKilled = sum(JF/epsilonF*B/m), 
      totEaten = sum(mortpred*B), 
      totGrowth=sum(Jtot*B/m)))  
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
