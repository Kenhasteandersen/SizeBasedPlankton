require(latex2exp)
source("basetools.R")

# ===================================================
# Plots for documentation
# ===================================================

plotAll = function() {
  pdfplot("AL.pdf", plotAL, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("AF.pdf", plotAF, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("Rstar.pdf", plotRstar, width=doublewidth, height = height)
}

calcC = function(volume, taxon) {
  C = 1e-6 * exp(-0.665) * volume^0.939 # Menden-deyer and Lessard (2000), table 4, for protists. mugC
  ixSmall = (volume<3000) & (!is.na(volume))
  C[ixSmall] = 1e-6 * exp(-0.583) * volume[ixSmall]^0.860
  ixDiatom = (taxon=="diatom") & (volume>3000) & (!is.na(volume))
  C[ixDiatom] = 1e-6 * exp(-0.933) * volume[ixDiatom]^0.881 # As above for diatom>3000 mum^3
  ixDiatom = (taxon=="diatom") & (volume<=3000)  & (!is.na(volume))
  C[ixDiatom] = 1e-6 * exp(-0.541) * volume[ixDiatom]^0.811 # As above for diatom>3000 mum^3
  return(C)
}

plotAL = function() {
  
  data = data.frame(C=NA,taxon=NA,AL=NA)
  #
  # Taguchi
  #
  #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
  #C = (Ata$C..pgC)*1e-6 # mugC
  #chl = 1e-3 * C/Ata$CperChl # mg chl
  #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
  #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
  #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
  #cat("AL = ", AL, "mugC/d/(wm2)\n")
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")

  C = calcC(Aed$volume, Aed$taxon)
  
  A = (Aed$alpha * 4.15) * C * (Aed$daylength/24)
  data = data.frame(C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
               
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  form = formula(log(A) ~ log( Cmax*C^(2/3) * a*C / (Cmax*C^(2/3) + a*C)))

  
  fit = nls(form,
            data = data,
            start = list(Cmax = .001, a=0.001),
            lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  fit_diatoms = nls(form,
            data = data[ixDiatom,],
            start = list(Cmax = .001, a=0.001),
            lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  fit_no_diatoms = nls(form,
            data = data[!ixDiatom,],
            start = list(Cmax = .001, a=0.001),
            lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  
  print(summary(fit_no_diatoms))
  #
  # Fits to 2/3 scaling
  #
  AL = exp(mean(log(data$A/data$C^(2/3))))
  AL_diatoms = exp(mean(log(data$A/data$C^(2/3))[ixDiatom]))
  AL_no_diatoms = exp(mean(log(data$A/data$C^(2/3))[!ixDiatom]))
  cat(AL_no_diatoms)
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim = data$C, ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
              xlab="Cell weight ($\\mu$gC)",
              ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/d/(W m$^2$)")
  points(data$C[!ixDiatom], data$A[!ixDiatom],pch=16, col="blue")
  points(data$C[ixDiatom], data$A[ixDiatom],pch=16, col="red")
  points(data$C[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="red")
  
  C = 10^seq(log10(min(data$C)), log10(max(data$C)), length.out = 100)
  lines(C, exp(predict(fit, list(C=C))), lwd=2)
  lines(C, exp(predict(fit_diatoms, list(C=C))), lwd=2, col="red")
  lines(C, exp(predict(fit_no_diatoms, list(C=C))), lwd=3, col="blue")

  lines(C, AL * C^(2/3), lty=dotted)
  lines(C, AL_diatoms * C^(2/3), col="red", lty=dotted)
  lines(C, AL_no_diatoms * C^(2/3), col="blue", lty=dotted)

  legend(x="bottomright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
         pch=16, col=c("red","Blue"))
}

plotALspecific = function() {
  data = data.frame(C=NA,taxon=NA,AL=NA)
  #
  # Taguchi
  #
  #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
  #C = (Ata$C..pgC)*1e-6 # mugC
  #chl = 1e-3 * C/Ata$CperChl # mg chl
  #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
  #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
  #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
  #cat("AL = ", AL, "mugC/d/(wm2)\n")
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  
  C = calcC(Aed$volume, Aed$taxon)
  
  A = (Aed$alpha * 4.15) * (Aed$daylength/24)
  data = data.frame(C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
  
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  form = formula(log(A) ~ log( Cmax*C^(-1/3) * a / (Cmax*C^(-1/3) + a)))
  
  fit = nls(form,
            data = data,
            start = list(Cmax = .001, a=0.001),
            lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  fit_diatoms = nls(form,
                    data = data[ixDiatom,],
                    start = list(Cmax = .001, a=0.001),
                    lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  fit_no_diatoms = nls(form,
                       data = data[!ixDiatom,],
                       start = list(Cmax = .001, a=0.001),
                       lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  
  print(summary(fit_no_diatoms))
  #
  # Fits to 2/3 scaling
  #
  AL = exp(mean(log(data$A/data$C^(-1/3))))
  AL_diatoms = exp(mean(log(data$A/data$C^(-1/3))[ixDiatom]))
  AL_no_diatoms = exp(mean(log(data$A/data$C^(-1/3))[!ixDiatom]))
  cat(AL_no_diatoms)
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim = data$C, ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
              xlab="Cell weight ($\\mu$gC)",
              ylab="Specific affinity for light, $\\textit{A_L}$ (1/d/(W m$^2$)")
  points(data$C[!ixDiatom], data$A[!ixDiatom],pch=16, col="blue")
  points(data$C[ixDiatom], data$A[ixDiatom],pch=16, col="red")
  points(data$C[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="red")
  
  C = 10^seq(log10(min(data$C)), log10(max(data$C)), length.out = 100)
  lines(C, exp(predict(fit, list(C=C))), lwd=2)
  lines(C, exp(predict(fit_diatoms, list(C=C))), lwd=2, col="red")
  lines(C, exp(predict(fit_no_diatoms, list(C=C))), lwd=3, col="blue")
  
  lines(C, AL * C^(-1/3), lty=dotted)
  lines(C, AL_diatoms * C^(-1/3), col="red", lty=dotted)
  lines(C, AL_no_diatoms * C^(-1/3), col="blue", lty=dotted)
  
  legend(x="topright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
         pch=16, col=c("red","Blue"))
}

plotALvolume = function() {
  data = data.frame(r=NA, V=NA,taxon=NA,AL=NA)
  #
  # Taguchi
  #
  #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
  #C = (Ata$C..pgC)*1e-6 # mugC
  #chl = 1e-3 * C/Ata$CperChl # mg chl
  #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
  #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
  #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
  #cat("AL = ", AL, "mugC/d/(wm2)\n")
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  C = calcC(Aed$volume, Aed$taxon)
  
  A = (Aed$alpha * 4.15) * (Aed$daylength/24)
  data = data.frame(V=Aed$volume, r = (3/(4*pi)*Aed$volume)^(1/3), 
                    Aed$taxon, A=A, source="Edwards", taxon=Aed$taxon,
                    mumax=Aed$mu_max)#rbind(data, 
  
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$r)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  form = formula(log(A) ~ log( cL*r^-1*(1-exp(-theta*r))))
  
  fit = nls(form,
            data = data,
            start = list(cL =1e-1, theta=1),
            lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
  fit_diatoms = nls(form,
                    data = data[ixDiatom,],
                    start = list(cL =1e-1, theta=1),
                    lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
  fit_no_diatoms = nls(form,
                       data = data[!ixDiatom,],
                       start = list(cL =1e-1, theta=1),
                       lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
  
  print(summary(fit_no_diatoms))
  #
  # Fits to -1 scaling
  #
  AL = exp(mean(log(data$A/data$r^(-1))))
  AL_diatoms = exp(mean(log(data$A/data$r^(-1))[ixDiatom]))
  AL_no_diatoms = exp(mean(log(data$A/data$r^(-1))[!ixDiatom]))
  cat(AL_no_diatoms)
  #
  # Plot
  #
  defaultplotvertical(2)
  loglogpanel(xlim = data$r, ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
              xaxis=FALSE,
              ylab="Affinity for light, $\\textit{A_L}$ ((d W m$^2$)$^{-1}$")
  points(data$r[!ixDiatom], data$A[!ixDiatom],pch=16, col="blue")
  points(data$r[ixDiatom], data$A[ixDiatom],pch=16, col="red")
  
  r = 10^seq(log10(min(data$r)), log10(max(data$r)), length.out = 100)
  lines(r, exp(predict(fit, list(r=r))), lwd=2)
  lines(r, exp(predict(fit_diatoms, list(r=r))), lwd=2, col="red")
  lines(r, exp(predict(fit_no_diatoms, list(r=r))), lwd=3, col="blue")
  
  lines(r, AL * r^(-1), lty=dotted)
  lines(r, AL_diatoms * r^(-1), col="red", lty=dotted)
  lines(r, AL_no_diatoms * r^(-1), col="blue", lty=dotted)
  
  legend(x="topright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
         pch=16, col=c("red","Blue"))
  #
  # Imax
  #
  fit_mu_diatoms = lm(log(mumax) ~ log(r), data=data[ixDiatom,])
  fit_mu_no_diatoms = lm(log(mumax) ~ log(r), data=data[!ixDiatom,])
  fit_mu = lm(log(mumax) ~ log(r), data=data)
  
  loglogpanel(xlim = data$r, ylim = data$mumax,
              xlab="radius ($\\mu$m)",
              ylab="max growth rate ($d^{-1}$)")
  points(data$r[!ixDiatom], data$mumax[!ixDiatom], pch=16, col="blue")
  points(data$r[ixDiatom], data$mumax[ixDiatom], pch=16, col="red")
  
  lines(r, exp(predict(fit_mu_diatoms, list(r=r))), lwd=2, col="red")
  lines(r, exp(predict(fit_mu, list(r=r))), lwd=2)
  lines(r, exp(predict(fit_mu_no_diatoms, list(r=r))), lwd=2, col="blue")
  #
  # Imax and alpha
  #
  defaultplot()
  x = data$A/exp(predict(fit, list(r=data$r)))
  y = data$mumax/exp(predict(fit_mu, list(r=data$r)))
  loglogpanel(xlim=x, ylim=x,
              xlab="alpha corrected", ylab="mumax corrected")
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  lines(c(0.01,100), c(0.01,100),lty=dotted)       
}

plotMumax = function() {
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  Aed$C = calcC(Aed$volume, Aed$taxon)
  #
  # Sort our nans:
  #
  Aed = Aed[!is.na(Aed$C), ]
  ixDiatom = Aed$taxon=="diatom"
  #
  # Fits
  #
  fit_diatoms = lm(log(mu_max) ~ log(C), data=Aed[ixDiatom,])
  fit_no_diatoms = lm(log(mu_max) ~ log(C), data=Aed[!ixDiatom,])
    
  print(fit_diatoms)
  print(fit_no_diatoms)
  
  defaultplot()
  loglogpanel(xlim = Aed$C, ylim = Aed$mu_max,
              xlab="Cell weight ($\\mu$gC)",
              ylab="max growth rate ($d^{-1}$)")
  points(Aed$C[!ixDiatom], Aed$mu_max[!ixDiatom], pch=16, col="blue")
  points(Aed$C[ixDiatom], Aed$mu_max[ixDiatom], pch=16, col="red")
  
  C = 10^seq(log10(min(Aed$C)), log10(max(Aed$C)), length.out = 100)
  lines(C, exp(predict(fit_diatoms, list(C=C))), lwd=2, col="red")
  lines(C, exp(predict(fit_no_diatoms, list(C=C))), lwd=2, col="blue")
}

plotMuAlphaCorrelation = function() {
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  Aed$C = calcC(Aed$volume, Aed$taxon)
  
  A = (Aed$alpha * 4.15) * C * (Aed$daylength/24)
  data = data.frame(C=C, taxon=Aed$taxon, A=A, mumax=Aed$mu_max)
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  ixDiatom = data$taxon=="diatom"
  data = data[ixDiatom,]
  #
  # Fits:
  #
  fit_mu = lm(log(mumax) ~ log(C), data=data)
  form = formula(log(A) ~ log( Cmax*C * a*C^(2/3) / (Cmax*C + a*C^(2/3))))
  fit_alpha = nls(form,
            data = data,
            start = list(Cmax = .001, a=0.001),
            lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  #
  # Plot corrected values
  #
  defaultplot()
  loglogpanel(xlim=c(0.1,10),ylim=c(0.1,10),
              xlab="Corrected alpha",
              ylab="Corrected mu_max")
  x = data$A/exp(predict(fit_alpha, list(C=data$C)))
  y = data$mumax/exp(predict(fit_mu, list(C=data$C)))
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  lines(c(0.1,10), c(0.1,10),lty=dotted)       
}

plotAF = function() {
  dat <- read.csv("../data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
  data = data.frame(w=1e3*dat$Body.mass, beta=24*0.001*dat$Fmax.1, Group=dat$Group)  
  
  ixProtist = (data$Group=="Nanoflagellates") | 
    (data$Group=="Dinoflagellates") | 
    (data$Group=="Dinoflagellate") | 
    (data$Group=="Ciliates") | 
    (data$Group=="Ciliate")
    
  x = data$beta/data$w
  x = x[ixProtist]
  x = x[!is.na(x)]
  SpecificBeta = exp(mean(log(x)))
  cat(SpecificBeta)
  
  defaultplot()
  loglogpanel(xlim=c(1e-6, 1), ylim=c(1e-8,1e-2),
              xlab = "Mass ($\\mu$gC)", ylab="Clearance rate $\\textit{A_F}$ (l/d)")
  points(data$w[ixProtist], data$beta[ixProtist])
  lines(data$w[ixProtist], SpecificBeta*data$w[ixProtist], lwd=2)
}

plotAN = function() {
  dat <- read.csv("../data/Nutrient traits from Edwards et al (2015b).csv", 
                  header=TRUE, sep=",")
  AN <- dat$vmax_nit / dat$k_nit
  m = calcC(dat$volume, dat$taxon)
  
  ix = !is.na(dat$c_per_cell)
  m[ix] = dat$c_per_cell[ix] # Use primary measurement if available
  
  defaultplot()
  loglogpanel(xlim=m, ylim=AN)
  points(m, AN)
}

plotRstar = function() {
  p = parameters()
  
  m = 10^seq(log10(0.3*min(p$m)), log10(max(p$m)), length.out = 1000)
  nu = pmin(1,p$c * m^(-1/3))
  mu = 0.4#c(0.2, 0.4, 0.6)
  
  tightaxes()
  defaultplot()
  loglogpanel(xlim=m, ylim=c(0.1,1000),
              xlab = "Cell mass ($\\mu$gC)",
              ylab = "Limiting concentration ($\\mu$g/l or W/m$^2$)")
  
  for (i in 1:length(mu)) {
    Nstar = m^0.667 * p$alphaJ/(p$AN*p$rhoCN) * mu[i]*(1-nu)/(p$alphaJ*(1-nu)-mu[i])
    Nstar[Nstar<0] = Inf
    lines(m, Nstar, lwd=2, col="blue")
    
    DOCstar = m^0.667 * p$alphaJ/(p$AN) * mu[i]*(1-nu)/(p$alphaJ*(1-nu)-mu[i])
    DOCstar[DOCstar<0] = Inf
    lines(m, DOCstar, lwd=2, col="magenta")
    
    #Lstar = m^0.333 * p$alphaJ/(p$AL) * mu[i]/(p$alphaJ*(1-nu)-mu[i])
    Lstar = p$alphaJ*(p$AL+p$cL*m^(1/3))*mu[i] / (p$AL*p$cL*(p$alphaJ-mu[i])*(1-nu))
    Lstar[Lstar<0] = Inf
    lines(m , Lstar, lwd=2, col="green")
  }
  
  legend(x="bottomright", bty="n",
         legend=c(TeX("$N^*$"), TeX("$DOC^*$"), TeX("$L^*$")),
         col=c("magenta","blue","green"),
         lwd=2)
}

calcQuadraticMortality <- function() {
  p = parameters()
  f0 = 0.6
  delta = p$m[2]/p$m[1]
  (1-f0)/delta * p$AF * sqrt(2*pi) * p$beta*p$sigma * exp(p$sigma^2/2)
}

