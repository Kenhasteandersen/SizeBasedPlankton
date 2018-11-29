library(deSolve)


#--------------------------------------------------
# Core logic for the model
#--------------------------------------------------

parameters <- function() {
  p = list()
  
  p$n = 10; # No of groups
  p$m = 10^seq(-8,1,length.out = p$n)  # Mass bins in mugC
  
  p$rhoCN = 5.68 # C:N ratio
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
  p$Jmax = 1 * p$m * (1-nu) # mugC/day
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
