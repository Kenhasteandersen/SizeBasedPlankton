require(latex2exp)
require(lattice)
source("basetools.R")
source("model.R")
source("modelChemostat.R")

# ===================================================
# Plots for documentation
# ===================================================

plotAll = function() {
  pdfplot("../aL.pdf", plot_aL, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../AF.pdf", plotAF, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../aN.pdf", plot_aN, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Mumax.pdf", plotMumax, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Rstar.pdf", plotRstar, width=doublewidth, height = height)
  pdfplot("../Mumax_corrected.pdf", plotMuAlphaCorrelation, width = 1.5*singlewidth, height=1.5*height)
  
  plotSimulationExamples()
  
  fontsize = trellis.par.get("fontsize")
  fontsize$text = 10
  trellis.par.set("fontsize",fontsize)
  plt=plotVaryLightAndDiffusion(n=15)
  pdf(file="../VaryLightAndDiffusion.pdf", width=doublewidth+0.5, height=height+.3)
  plt
  dev.off()
  
  pdfplot("../Functions.pdf", plotFunctions,n=25, width=doublewidth, height=2*height)
  
  pdfplot("../GridPreference.pdf", plotGridPreference, width=1.5*singlewidth, height=height)
  
  pdfplot("../Gridtest.pdf", plotGridtest, width=doublewidth, height=height)
  
  system("cp ../*pdf ../dropbox")
  #system2("cp *pdf ~/Dropbox/Apps/Overleaf/A\ minimal\ size-based\ model\ of\ unicellular\ plankton/")
}

convertVolume2Mass = function(vol, taxon="other") {
  C = 1e-6 * exp(-0.665) * vol^0.939 # Menden-deyer and Lessard (2000), table 4, for protists. mugC
  ixSmall = (vol<3000) & (!is.na(vol))
  C[ixSmall] = 1e-6 * exp(-0.583) * vol[ixSmall]^0.860
  ixDiatom = (taxon=="diatom") & (vol>3000) & (!is.na(vol))
  C[ixDiatom] = 1e-6 * exp(-0.933) * vol[ixDiatom]^0.881 # As above for diatom>3000 mum^3
  ixDiatom = (taxon=="diatom") & (vol<=3000)  & (!is.na(vol))
  C[ixDiatom] = 1e-6 * exp(-0.541) * vol[ixDiatom]^0.811 # As above for diatom>3000 mum^3
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
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  
  A = Aed$alpha * C * (Aed$daylength/24) # convert to units of mu gC/(mu mol photons/m2/s),
  # corrected for daylength.
  data = data.frame(C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
  
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
  form = formula(log(A) ~ log( a*C^(2/3) * (1 - exp(-Cmax*C^(1/3)) ) ))
  
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
  loglogpanel(xlim = data$C, 
              ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
              xlab="Cell weight ($\\mu$gC)",
              ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/($\\mu$ mol photons m$^{-2}s^{-1})")
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
  
  # Camila's parameters:
  cL=0.08733
  AL=0.004864  
  lines(C, cL*C * AL*C^(2/3) / ( cL*C + AL*C^(2/3) ) , col="orange")
  
  legend(x="bottomright", bty="n",
         legend=c("Diatoms", "Other phototrophs", "Camila"),
         pch=c(16,16,NA),
         lwd=c(0,0,2),
         col=c("red","Blue","Orange"))
}

plotALsimple = function() {
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
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  r = (Aed$volume*3/(4*pi))^(1/3)
  
  A = Aed$alpha * C * (Aed$daylength/24) # convert to units of mu gC/day/(mu mol photons/m2/s),
  # corrected for daylength.
  data = data.frame(r=r, C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
  
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
  form = formula(log(A) ~ log( AL*r^2 * (1 - exp(-cL*C^(1/3)) ) ))
  
  fit = nls(form,
            data = data,
            start = list(cL = 1, AL=1e-7),
            lower=list(cL=1e-20, AL=1e-20), algorithm="port")

  print(summary(fit))
  cat(c("r* = ",1/(coef(fit)[[1]]*(0.19e-6*4/3*pi)^(1/3)) ), "mu m\n")
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
  loglogpanel(xlim = c(0.1, 300), 
              ylim = c(0.5*min(data$A[!is.na(data$A)]), 2*max(data$A[!is.na(data$A)])),
              xlab="Cell radius ($\\mu$m)",
              ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/day/($\\mu$ mol photons m$^{-2}s^{-1})")
  points(data$r[!ixDiatom], data$A[!ixDiatom],pch=15, col="darkgreen")
  points(data$r[ixDiatom], data$A[ixDiatom],pch=16, col="darkgreen")
  points(data$r[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="darkgreen")
  
  r = 10^seq(-1, 4, length.out = 100)
  C = 4/3*pi*r^3*0.2e-6
  lines(r, exp(predict(fit, list(r=r,C=C))), lwd=2)

  lines(r, coef(fit)[[2]]*r^2, lty=dotted)
  lines(r, coef(fit)[[2]]*r^2*coef(fit)[[1]]*C^(1/3), lty=dotted)

  legend(x="bottomright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
           pch=c(15,16),
         lwd=0,
         col="darkgreen")
}

plot_aL = function() {
  data = data.frame(r=NA,C=NA,taxon=NA,aL=NA)
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
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  r = (Aed$volume*3/(4*pi))^(1/3)
  
  aL = Aed$alpha * (Aed$daylength/24) # convert to units of mu gC/day/(mu mol photons/m2/s),
  # corrected for daylength.
  data = data.frame(r=r, C=C, taxon=Aed$taxon, aL=aL, source="Edwards")#rbind(data, 
  
  #
  # Sort out nans:
  #
  data = data[(!is.na(data$aL) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
  form = formula(log(aL) ~ log( alphaL/r * (1 - exp(-r/rStar) ) ))
  
  fit = nls(form,
            data = data,
            start = list(rStar = 1.5, alphaL=1),
            lower=list(rStar=0, alphaL=1e-20), 
            upper=list(rStar=10, alphaL=1e10), algorithm="port")
  
  print(summary(fit))
  cat(c("alpha_L = ", coef(fit)[[2]], "1/day/mugC*mum/(light)"))
  cat(c("r* = ", coef(fit)[1], "mu m\n"))
  #
  # Fits to 2/3 scaling
  #
  #alphaL = exp(mean(log(data$A/data$C^(2/3))))
  #AL_diatoms = exp(mean(log(data$A/data$C^(2/3))[ixDiatom]))
  #AL_no_diatoms = exp(mean(log(data$A/data$C^(2/3))[!ixDiatom]))
  #cat(AL_no_diatoms)
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim = c(0.1, 300), 
              ylim = c(0.5*min(data$aL[!is.na(data$aL)]), 2*max(data$aL[!is.na(data$aL)])),
              xlab="Cell radius ($\\mu$m)",
              ylab="Affinity for light, $\\textit{a_L}$ (day$\\cdot \\mu$ mol photons m$^{-2}s^{-1}\\cdot \\mu$gC)$^{-1}$")
  points(data$r[!ixDiatom], data$aL[!ixDiatom],pch=15, col="darkgreen")
  points(data$r[ixDiatom], data$aL[ixDiatom],pch=16, col="darkgreen")
  points(data$r[data$source=="Taguchi"], data$aL[data$source=="Taguchi"], pch=17, col="darkgreen")
  
  r = 10^seq(-1, 4, length.out = 100)
  lines(r, exp(predict(fit, list(r=r,C=C))), lwd=2)
  
  lines(r, coef(fit)[[2]]/r, lty=dotted)
  lines(r, coef(fit)[[2]]/coef(fit)[[1]]*r/r, lty=dotted)
  
  legend(x="topright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
         pch=c(15,16),
         lwd=0,
         col="darkgreen")
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
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  
  A = (Aed$alpha) * (Aed$daylength/24)
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
              ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/day/($\\mu$ mol photons m$^{-2}s^{-1})")
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
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  Aed$C =   C
  # Sort our nans:
  Aed = Aed[!is.na(Aed$C), ]
  # Correct to 10 degrees with a Q10 of 1.5
  Aed$mu_max = Aed$mu_max * 1.5^((10-Aed$temperature)/10)
  
  A = data.frame(species=Aed$species, taxon=as.factor(Aed$taxon), 
                 C=Aed$C, mu_max=Aed$mu_max)
  #
  # Kiørboe and hirst (2014):
  #
  Akh = read.csv('../Data/Kiorboe and Hirst (2014) maximum growth rates.csv', 
                 sep=",", header=TRUE, skip=1, as.is=TRUE)
  Akh$Group[Akh$Group=="Dinoflagellates"] = "dinoflagellate"
  Akh$Group[Akh$Group=="Ciliates"] = "ciliate"
  # Correct from 15 to 10 degrees with a Q10 of 2.8
  Akh$Specific.growth.1 = Akh$Specific.growth.1 * 2.8^((10-15)/10)
  
  A = rbind(A, 
            data.frame(species=Akh$Species, taxon=as.factor(Akh$Group),
                       C=Akh$Body.Mass*1000, mu_max=Akh$Specific.growth.1*24))
  #
  # Kirchman (2016) - but not max growth rates
  #
  Ak = read.csv('../Data/ma08_kirchman_supappendix3.csv',
                sep=",", header=TRUE, skip=15)
  Ak = Ak[1:149,]
  # Sort out incomplete entries:
  Ak = Ak[Ak[,7]!="Used biovol*",]
  Ak = Ak[!is.na(Ak[,5]),]

  A = rbind(A,
            data.frame(species="NA",
                       taxon="bacteria", 
                       C=as.numeric(as.character(Ak$fgC.cell))*1e-9, 
                       mu_max=as.numeric(Ak$Rate...d.)))
  #
  # Rose and Caron (2007)
  #
  Arc = read.csv('../Data/Rose and Caron bacterivores.csv',
                 sep=",", header=TRUE)
  Arc$volume = as.numeric( gsub(",","", as.character(Arc$volume)) )
  # Correct to 10 degrees with a Q10 of 2.8
  Arc$Growth.Rate = Arc$Growth.Rate * 2.8^((10-Arc$Temperature)/10)
  
  A = rbind(A, data.frame(
    species="NA",
    taxon="bacterivore",
    C = convertVolume2Mass(Arc$volume),
    mu_max=as.numeric(Arc$Growth.Rate)))
  
  Arc = read.csv('../Data/Rose and Caron herbivores.csv',
                 sep=",", header=TRUE)
  Arc$volume = as.numeric( gsub(",","", as.character(Arc$volume)) )
  # Correct to 10 degrees with a Q10 of 2.8
  Arc$Growth.Rate = Arc$Growth.Rate * 2.8^((10-Arc$Temperature)/10)
  
  A = rbind(A,data.frame(
    species="NA",
    taxon="herbivore",
    C = convertVolume2Mass(Arc$volume),
    mu_max = as.numeric(Arc$Growth.Rate)))
  #
  # Define groups:
  #
  ixDiatom = A$taxon=="diatom"
  ixMixotroph = A$taxon=="nanoflagellate" | A$taxon=="dinoflagellate" 
  ixHeterotroph = A$taxon=="ciliate" 
  ixBacteria = A$taxon=="bacteria"
  ixPhototroph = !ixDiatom & !ixMixotroph & !ixHeterotroph &!ixBacteria &
    A$taxon!="bacterivore" & A$taxon!="herbivore"
  #
  # Fits
  #
  fit_diatoms = lm(log(mu_max) ~ log(C), data=A[ixDiatom,])
  fit_mixotrophs = lm(log(mu_max) ~ log(C), data=A[ixMixotroph,])
  fit_heterotrophs = lm(log(mu_max) ~ log(C), data=A[ixHeterotroph,])
  fit_phototrophs = lm(log(mu_max) ~ log(C), data=A[ixPhototroph,])
  
  fit = nls(mu_max ~ (1-parameters()$c * C^(-1/3))*mu/(mu*c*C^0.3333+1), 
            data=A[!ixBacteria,],
            start=list(mu=1, c=1))
  
  defaultplot()
  semilogxpanel(xlim = c(1e-9,1), ylim = A$mu_max,
              xlab="Cell weight ($\\mu$gC)",
              ylab="max growth rate ($d^{-1}$)")
  points(A$C[ixMixotroph], A$mu_max[ixMixotroph], pch=16, col="blue")
  points(A$C[A$taxon=="bacterivore"], A$mu_max[A$taxon=="bacterivore"], pch=1, col="blue", cex=.1)
  points(A$C[ixHeterotroph], A$mu_max[ixHeterotroph], pch=16, col="red")
  points(A$C[A$taxon=="herbivore"], A$mu_max[A$taxon=="herbivore"], pch=1, col="red", cex=0.1)
  points(A$C[ixPhototroph], A$mu_max[ixPhototroph], pch=16, col="green")
  points(A$C[ixDiatom], A$mu_max[ixDiatom], pch=16, col="darkgreen")
  points(A$C[ixBacteria], A$mu_max[ixBacteria], pch=1, col="brown", cex=0.1)
  
  
  C = 10^seq(-9, log10(max(A$C)), length.out = 100)
  #lines(C, exp(predict(fit_diatoms, list(C=C))), col="darkgreen")
  #lines(C, exp(predict(fit_mixotrophs, list(C=C))), col="blue")
  #lines(C, exp(predict(fit_heterotrophs, list(C=C))), col="red")
  #lines(C, exp(predict(fit_phototrophs, list(C=C))), col="green")
  
  m = 10^seq(-9,1,length.out = 100)
  lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)), lwd=3)
  #lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)) - parameters()$cR, lwd=2)
  #lines(C, predict(fit, list(C=C)))
  
  # Add the used curve:
  p = parameters()
  #lines(p$m, p$Jmax/p$m, col="black", lwd=2)
  
  # Maximum uptake of phagotrophy:
  lines(m, p$epsilonF*p$cF*m^(-1/3)-p$cR*p$alphaJ, col="red", lwd=2)

  # Camila
  #lines(C, 0.12*C^-0.25, col="orange", lwd=1)
  #lines(C, (0.12-0.03)*C^-0.25, col="orange", lwd=2)
  
  legend(x="topleft", bty="n",
                  legend=c("Diatoms","Other phototrophs",
                           "Mixotrophs", "Heterotrophs",
                           "Bacteria","Model","Max. phagotrophy"),
                  lty=c(NA,NA,NA,NA,NA,solid,solid),
                  pch=c(16,16,16,16,16,NA,NA),
                  lwd=c(0,0,0,0,0,2,2),
                  col=c("darkgreen","green","blue","red","brown","black","red"))
}

plotMuAlphaCorrelation = function() {
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  Aed$C = C
  
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
  # Plot corrected values (using residuals)
  #
  defaultplot(mfcol=c(1,2))
  loglogpanel(xlim=c(0.1,10),ylim=c(0.1,10),
              xlab="Residual light-affinity",
              ylab="Residual $\\mu_{max}$")
  x = data$A/exp(predict(fit_alpha, list(C=data$C)))
  y = data$mumax/exp(predict(fit_mu, list(C=data$C)))
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  #lines(c(0.1,10), c(0.1,10),lty=dotted)     

  fit = lm(log(x) ~ log(y))
  xx = 10^seq(-1,2,length.out = 100)
  #lines(xx, exp(fit$coefficients[1])*xx^fit$coefficients[2])
  legend(x="topleft", pch=16, col=c("red","blue"),
         legend=c("Diatoms","Others"))
  #
  # Plot mass-specific values
  #
  loglogpanel(xlim=c(0.001,2),ylim=c(0.1,5),
              xlab="Specific light-affinity",
              ylab="Specific $\\mu_{max}$")
  x = data$A/data$C
  y = data$mumax
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  #lines(c(0.1,10), c(0.1,10),lty=dotted)     
  legend(x="topleft", pch=16, col=c("red","blue"),
         legend=c("Diatoms","Others"))
  
  fit = lm(log(x) ~ log(y))
  xx = 10^seq(-3,1,length.out = 100)
  #lines(xx, exp(fit$coefficients[1])*xx^fit$coefficients[2])
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
  
  # Camila
  C = range(data$w[ixProtist])
  lines(C, 0.018*C, col="orange")
  
  
}
#
# Plot specific affinity:
#
plot_aN = function() {
  dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
                 header=TRUE,sep=",")
  #
  # Convert carbon from mol to g:
  #
  dat$c_per_cell = dat$c_per_cell*12
  #
  # Convert to weights when only volume is given:
  # 
  ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
  dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
  C = dat$c_per_cell
  #
  # Convert to ESD/2 (radius) assuming spherical cells
  #
  r = ( 3/(4*pi)*dat$volume )^(1/3)
  
  #
  # Calc affinities
  #
  a_nit = dat$vmax_nit / dat$k_nit / dat$c_per_cell
  a_amm = dat$vmax_amm / dat$k_amm / dat$c_per_cell 
  a_phos = dat$vmax_p / dat$k_p / dat$c_per_cell
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim=c(0.1 ,200), ylim=c(1e-5, 20), #c(1e-9 ,1)
              xlab="Radius ($\\mathrm{\\mu}$m)",
              ylab="Nutrient affinity, $\\textit{a}_N$ (L/day/$\\mathrm{\\mu}$gC)")
  
  col = 1
  taxons = unique(dat$taxon[!is.na(C) & !(is.na(Anit) & is.na(Aamm)  & is.na(Aphos))])
  aff = data.frame(r=NULL, A=NULL)
  for (i in taxons) {
    ix = dat$taxon == i
    points(r[ix], a_nit[ix], pch=16, col=col)
    points(r[ix], a_amm[ix], pch=17, col=col)
    points(r[ix], a_phos[ix], pch=18, col=col)
    
    #points(C[ix], a_nit[ix], pch=16, col=col)
    #points(C[ix], a_amm[ix], pch=17, col=col)
    #points(C[ix], a_phos[ix], pch=18, col=col)
    col = col + 1
    
    aff = rbind( aff, data.frame(r=r[ix], a= a_nit[ix]))
    aff = rbind( aff, data.frame(r=r[ix], a= a_amm[ix]))
    aff = rbind( aff, data.frame(r=r[ix], a= a_phos[ix]))
  }
  ix = !is.na(aff$r) & !is.na(aff$a)
  aff = aff[ix,]
  #points(C, Aphos,pch=18)
  legend(x="topright",bty="n",
         legend=taxons, pch=16, col=seq(1,length(taxons)))
  
  r = 10^seq(-1,3,length.out = 100) # mum
  Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
  rho = 0.57 # g/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )
  aNmax = 3*Diff*1e8*1e-3/rho * 1e-6 # L/day/mugC
  rstar = 2; # mum
  corr = 1 - parameters()$c*m^(-1/3)
  #cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
  #lines(m, ANmax/m, lwd=1, lty=dotted)  
  #lines(m, 0.2*m/m, lwd=1, lty=dotted)  
  lines(r, aNmax*r^-2, lwd=1, lty=dotted)  
  #lines(m, 0.2*m/m, lwd=1, lty=dotted)  
  lines(r, aNmax*r^-2*(r/rstar)^2, lwd=1, lty=dotted, col="blue")
  lines(r, aNmax*r^-2/(1+((r)/rstar)^-2), lwd=2)
  
  #lines(m, p$AN/p$cN*m^(2/3)/m, lty=dotted, col="blue")
  #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
  cat("aN = ", aNmax, "r^-2 \n")
  cat("r_D^* = ", rstar," mum\n")
  
  # Make a bi-linear fit:
  #form = formula( log(C) ~ log(a/(1+b*C^-0.333))+0.3333*log(C) )
  #fit = nls(form, data=aff, 
  #          start = list(a=1e-4, b=0.05),
  #          lower = list(a=0,b=0), algorithm = "port", trace=TRUE)
  #lines(m, exp(predict(fit,list(C=m))))
  
  #V = 10^seq(-2,8,length.out = 100)
  #alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
  #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
  
  # Camila's parameters:
  #lines(C, 3.75e-5*C^(1/3), col="orange")
  
  #alphaN_Banas2011 = 2.6/0.1 / 14 / 6 # convert from uM N-> ugN, from gN->gC
  #lines(m , alphaN_Banas2011*m^(1-0.45), col="blue", lty=dashdotted)
  
  legend(x="bottomleft", bty="n",
         legend=c("Model","Diffusion limit","Porter limit"),
         lty=c(solid,dotted,dotted),
         lwd=c(2,1,1,1),
         col=c("black","black","blue","orange"))
  
  #r = (3*dat$volume/(4*pi))^(1/3) # mu m
  #m = 0.3e-6*(2*r)^3
  #ANmax = 4*pi*Diff*(r*1e-4) * 1e-3
  #points(m, ANmax, col="red")
}



plotAN = function() {
  dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
                 header=TRUE,sep=",")
  #
  # Convert carbon from mol to g:
  #
  dat$c_per_cell = dat$c_per_cell*12
  #
  # Convert to weights when only volume is given:
  # 
  ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
  dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
  C = dat$c_per_cell
  #
  # Calc affinities
  #
  Anit = dat$vmax_nit / dat$k_nit
  Aamm = dat$vmax_amm / dat$k_amm
  Aphos = dat$vmax_p / dat$k_p
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim=C, ylim=c(1e-9,1e-3),
              xlab="Mass ($\\mathrm{\\mu g_C}$)",
              ylab="Nutrient affinity, $\\textit{A_N}$ (L/day)")
  
  col = 1
  taxons = unique(dat$taxon[!is.na(C) & !(is.na(Anit) & is.na(Aamm)  & is.na(Aphos))])
  aff = data.frame(C=NULL, A=NULL)
  for (i in taxons) {
    ix = dat$taxon == i
    points(C[ix], Anit[ix], pch=16, col=col)
    points(C[ix], Aamm[ix], pch=17, col=col)
    points(C[ix], Aphos[ix], pch=18, col=col)
    col = col + 1
    
    aff = rbind( aff, data.frame(C=C[ix], A= Anit[ix]))
    aff = rbind( aff, data.frame(C=C[ix], A= Aamm[ix]))
    aff = rbind( aff, data.frame(C=C[ix], A= Aphos[ix]))
  }
  ix = !is.na(aff$C) & !is.na(aff$A)
  aff = aff[ix,]
  #points(C, Aphos,pch=18)
  legend(x="bottomright",bty="n", 
         legend=taxons, pch=16, col=seq(1,length(taxons)))
  
  
  m = 10^seq(-7,1) # mu gC
  r = 1.5/2*(1e-6*m)^(1/3) # cm
  Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
  ANmax = 4*pi*Diff*r*1e-3 # L/day
  cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
  lines(m, ANmax, lwd=1, lty=dotted)  
  lines(parameters()$m, parameters()$ANm, lwd=2)
  lines(m, p$AN/p$cN*m^(2/3), lty=dotted, col="blue")
  #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
  cat("alphaN = ", mean(parameters()$AN*m^(1/3)), "\n")
  
  # Make a bi-linear fit:
  #form = formula( log(C) ~ log(a/(1+b*C^-0.333))+0.3333*log(C) )
  #fit = nls(form, data=aff, 
  #          start = list(a=1e-4, b=0.05),
  #          lower = list(a=0,b=0), algorithm = "port", trace=TRUE)
  #lines(m, exp(predict(fit,list(C=m))))
  
  #V = 10^seq(-2,8,length.out = 100)
  #alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
  #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
  
  # Camila's parameters:
  #lines(C, 3.75e-5*C^(1/3), col="orange")
  
  #alphaN_Banas2011 = 2.6/0.1 / 14 / 6 # convert from uM N-> ugN, from gN->gC
  #lines(m , alphaN_Banas2011*m^(1-0.45), col="blue", lty=dashdotted)
  
  legend(x="topleft", bty="n",
         legend=c("Model","Diffusion limit","Porter limit"),
         lty=c(solid,dotted,dotted),
         lwd=c(2,1,1,1),
         col=c("black","black","blue","orange"))
  
  #r = (3*dat$volume/(4*pi))^(1/3) # mu m
  #m = 0.3e-6*(2*r)^3
  #ANmax = 4*pi*Diff*(r*1e-4) * 1e-3
  #points(m, ANmax, col="red")
}

plotANsimple = function() {
  dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
                 header=TRUE,sep=",")
  #
  # Convert carbon from mol to g:
  #
  dat$c_per_cell = dat$c_per_cell*12
  #
  # Convert to weights when only volume is given:
  # 
  ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
  dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
  C = dat$c_per_cell
  #
  # Calc affinities
  #
  Anit = dat$vmax_nit / dat$k_nit
  Aamm = dat$vmax_amm / dat$k_amm
  Aphos = dat$vmax_p / dat$k_p
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim=C, ylim=c(1e-9,1e-3),
              xlab="Mass ($\\mathrm{\\mu g_C}$)",
              ylab="Nutrient affinity, $\\textit{A_N}$ (L/day)")
  
  col = 1
  taxons = unique(dat$taxon[!is.na(C) & !(is.na(Anit) & is.na(Aamm)  & is.na(Aphos))])
  for (i in taxons) {
    ix = dat$taxon == i
    points(C[ix], Anit[ix], pch=16, col="darkblue")
    points(C[ix], Aamm[ix], pch=17, col="darkblue")
    #points(C[ix], Aphos[ix], pch=18, col=col)
    col = col + 1
  }
  #points(C, Aphos,pch=18)
  #legend(x="bottomright",bty="n", 
  #       legend=taxons, pch=16, col=seq(1,length(taxons)))
  
  
  m = 10^seq(-7,1) # mu gC
  r = 1.5/2*(1e-6*m)^(1/3) # cm
  Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
  ANmax = 4*pi*Diff*r*1e-3 # L/day
  cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
  lines(m, ANmax, lwd=1, lty=dotted)  
  lines(m, 1e-3*m^(2/3))
  #lines(m, parameters()$AN*m^(1/3), lwd=2)
  #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
  cat("alphaN = ", mean(parameters()$AN*m^(1/3)), "\n")
  
  V = 10^seq(-2,8,length.out = 100)
  alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
  #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
}

plotLimitation = function(p=parameters()) {
  m = p$m
  
  DOC = 2
  N = 1
  L = 60
  B = 10
  r = calcRates(L, N, DOC, rep(B, p$n), p)
  
  defaultplot()
  semilogxpanel(xlim=m, ylim=c(0, 1.5))
  
  lines(m, r$JN/m, col="blue")
  lines(m, r$JDOC/m, col="brown")
  lines(m, r$JL/m, col="darkgreen")
  lines(m, r$JF/m, col="red")
  lines(m, r$Jmax/m, col="black")
  lines(m, r$Jtot/m, col="black")
}

plotRstar = function() {
  p = parameters()
  
  m = 10^seq(log10(5e-9), log10(max(p$m)), length.out = 1000)
  nu = pmin(1,p$c * m^(-1/3))
  mu = c(0,0.4)#c(0.2, 0.4, 0.6)
  
  tightaxes()
  defaultplot()
  loglogpanel(xlim=m, 
              ylim=c(0.0005,1000),
              xlab = "Cell mass ($\\mu$gC)",
              ylab = "Limiting concentration ($\\mu$mol/l) or ($\\mu$E/m$^2$/s)")
  
  lty=dotted
  lwd=1
  for (i in 1:length(mu)) {
    #correct = mu[i]*(1-nu)/((1-nu)-mu[i])
    #Nstar = m^0.667 /(p$AN*p$rhoCN) * correct

    correct = p$cLeakage + p$alphaJ*mu[i]*(1-nu)/(p$alphaJ*(1-nu)-mu[i])
    Nstar = m^(2/3)/(p$AN*p$rhoCN) * correct
    Nstar[Nstar<0] = Inf
    lines(m, Nstar/14, lwd=lwd, lty=lty,col="blue")
    
    correct = p$cLeakage*m^(-1/3) + p$alphaJ*(p$cR+ mu[i]*(1-nu)/(p$alphaJ*(1-nu)-mu[i]))
    DOCstar = m^(2/3)/(p$AN) * correct
    DOCstar[DOCstar<0] = Inf
    lines(m, DOCstar/12, lwd=lwd, lty=lty,col="magenta")
    
    #Lstar = m^0.333 * p$alphaJ/(p$AL) * mu[i]/(p$alphaJ*(1-nu)-mu[i])
    #Lstar = p$alphaJ*(p$AL+p$cL*m^(1/3))*mu[i] / (p$AL*p$cL*(p$alphaJ-mu[i])*(1-nu))
    Lstar = m^(1/3)/(p$AL*p$epsilonL) * 1/(1 -exp(-p$cL*m^(1/3))) * correct
    Lstar[Lstar<0] = Inf
    lines(m , Lstar, lwd=lwd, lty=lty,col="darkgreen")
    
    
    Fstar = 1/(p$AF*p$epsilonF) * correct
    Fstar[Fstar<0] = Inf
    lines(m, Fstar/12, lwd=lwd, lty=lty, col="darkred")
    
    lty=1
    lwd=2
  }
  
  legend(x="bottomright", bty="n",
         legend=c(TeX("$\\textit{N}^*$"), 
                  TeX("$$DOC^*$"), 
                  TeX("$$\\textit{L}^*$"),
                  TeX("$$\\textit{F}^*$")),
         col=c("blue","magenta","darkgreen","darkred"),
         lwd=2)
}

calcMufactor = function(p) {
  f0 = 0.6
  delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
  cat("Delta = ", delta, "\n")
  return( (1-f0) * p$AF *sqrt(2*pi)*p$sigma )
}

plotSimulationExamples = function() {
  #
  # Oligotrophic
  #
  p = parametersChemostat()
  p$d = 0.0025
  p$L = 100
  p$mortHTL = 0.03
  
  sim = simulate(p)
  pdfplot("../spectrum1.pdf",FUN=plotSpectrum, sim=sim, width=doublewidth, height=1.5*height)
  pdfplot("../rates1.pdf", FUN=plotRates, sim=sim, width=doublewidth, height=1.5*height)
  #
  # Eutrophic
  #
  p = parametersChemostat()
  p$d = 0.3
  p$L = 100
  p$mortHTL = 0.03
  
  sim = simulate(p)
  pdfplot("../spectrum2.pdf",FUN=plotSpectrum, sim=sim, width=doublewidth, height=1.5*height)
  pdfplot("../rates2.pdf", FUN=plotRates, sim=sim, width=doublewidth, height=1.5*height)
}


plotVaryDiffusion = function(p=parametersChemostat(),
                             L=100,
                             d=10^seq(-6,-3,length.out=10)) {
  
  B = matrix(data=0, nrow=p$n, ncol=length(d))
  for (i in 1:length(d)) {
    p$d=d[i]
    sim = simulateChemostat(p)
    B[,i] = sim$B
  }
  image(log10(d),log10(p$m),t(B),
        xlab="log10(mixing)", ylab="log10(cell mass)")
}

plotVaryLightAndDiffusion = function(n=10) {
  library(lattice)
  library(latticeExtra)
  library(reshape)
  
  d = seq(0.004,.5,length.out=n) #10^seq(-2,log10(2),length.out = 10)
  L = seq(10,200,length.out=n)
  #
  # Simulations
  #
  df = data.frame(size=NULL, d=NULL, L=NULL, biomass=NULL)
  df2 = data.frame(func=NULL, d=NULL, L=NULL, z=NULL)
  #tmp = matrix(ncol=4, nrow=0)
  #B = array(dim=c(3,length(d), length(L)))
  for (i in 1:length(d))
    for (j in 1:length(L)) {
      p = parametersChemostat()
      p$d = d[i]
      p$L = L[j]
      print(p$d)
      print(p$L)
      #p$mHTL = 0.008
      sim = simulate(p)
      func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
      
      df = rbind(df, data.frame(size="Pico",d=d[i],L=L[j],biomass=func$Bpico))
      df = rbind(df, data.frame(size="Nano",d=d[i],L=L[j],biomass=func$Bnano))
      df = rbind(df, data.frame(size="Micro",d=d[i],L=L[j],biomass=func$Bmicro))
      #df = rbind(df, data.frame(size="Chl",d=d[i],L=L[j],biomass=func$Chl_per_l))
      
      #tmp = rbind(tmp, c(func$Bpico, func$Bnano, func$Bmicro, func$Chl_per_l))
      
      df2 = rbind(df2, data.frame(func="Biomass", d=d[i],L=L[j],
                                  z=func$Bpico+func$Bnano+func$Bmicro))
      df2 = rbind(df2, data.frame(func="N", d=d[i],L=L[j], z=sim$N))
      df2 = rbind(df2, data.frame(func="DOC", d=d[i],L=L[j], z=sim$DOC))
      
      #B[1,i,j] = func$Bpico+func$Bnano+func$Bmicro
      #B[2,i,j] = sim$N
      #B[3,i,j] = sim$DOC
    }
  #
  # Plots
  #
  plt=levelplot(biomass ~ d*L | size, data=df,
                aspect=1, region=TRUE, 
                col.regions = terrain.colors(100),
                #scales = list(x = list(log = 10, at=c(0.01,0.1,1))),
                #xscale.components = xscale.components.logpower,#10ticks(n=5,n2=10),
                layout=c(3,1),
                ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), xlab=TeX("Mixing rate (day$^{-1}$)"))
  
  plt2 = levelplot(z ~ d*L | func, data=df2,
                   aspect=1, region=TRUE,
                   col.regions = terrain.colors(100),
                   layout=c(3,1),
                   ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), 
                   xlab=TeX("Mixing rate (day$^{-1}$)"))
  
  return(plt)
}

plotFunctions = function(L=c(35, 100), n=10) {
  
  panelsFunctions = function(L=c(18,40), n=10, yaxis=TRUE) {
    d = 10^seq(-2,log10(.7),length.out = n) #seq(0.02,2,length.out=n) #
    
    p = parametersChemostat()
    p$L = L
    
    F = data.frame()
    for (i in 1:length(d)) {
      print(d[i])
      p$d = d[i]
      sim = simulate(p)
      func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
      func$d = d[i]
      func$N = sim$N
      func$DOC = sim$DOC
      
      F = rbind(F, as.data.frame(func))
    }
    
    # Biomass
    B = F$Bpico+F$Bnano+F$Bmicro
    defaultpanel(xlim=c(0, max(d)), ylim=c(0,0.15), xaxis=FALSE, yaxis=yaxis, bty="l",
                 ylab="Biomass (gC/m$^2$)")
    lines(F$d, B, lwd=2)
    lines(F$d, F$Bpico, lwd=0.5, col="grey")
    lines(F$d, F$Bnano, lwd=1, col="grey")
    lines(F$d, F$Bmicro, lwd=1.5, col="grey")
    
    
    # N
    semilogypanel(xlim=d, ylim=c(0.001,20), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="N ($\\mu$mol N/l)")
    lines(F$d, F$N/14, lwd=2)
    
    #DOC
    semilogypanel(xlim=d, ylim=c(1,5000), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="DOC  ($\\mu$mol C/l)")
    lines(F$d, 1000*F$DOC/12, lwd=2)
    
    #Production:
    semilogypanel(xlim=d, ylim=c(10,3000), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="Prod. (gC/m$^2$/yr)")
    lines(F$d, F$prodCgross)
    lines(F$d, F$prodCnet, col="blue")
    lines(F$d, F$prodNew, col="green")
    lines(F$d, F$prodHTL, col="red")
    
    # Eff
    defaultpanel(xlim=d, ylim=c(0,1), bty="l",
                 ylab="$\\epsilon_{HTL}$",
                 xlab="Mixing rate (day$^{-1}$)",
                 yaxis=yaxis)
    lines(F$d, F$effHTL, lwd=2)
  }
  
  defaultplotvertical(mfcol=c(5,2))
  
  panelsFunctions(L=L[1],n=20)
  panelsFunctions(L=L[2], yaxis=FALSE,n=20)
}


plotPreferences = function(p = parameters()) {
  defaultplot()
  semilogxpanel(xlim=p$m, xlab="Cell weight ($\\mu$gC)",
                ylim=c(0,1.1), ylab="Preference")
  
  m = 10^seq(-10,1,length.out = 1000)
  for (i in 1:p$n) {
    lines(m, phi(p$m[i]/m,p$beta,p$sigma), col=grey(0.8))
    lines(p$m, phi(p$m[i]/p$m,p$beta,p$sigma), lty=dotted)
    points(p$m, phi(p$m[i]/p$m,p$beta,p$sigma),pch=16)
  }
  
}

plotGridPreference = function(p = parameters()) {
  m = 10^seq(-6,1,length.out = 100)
  defaultplot()

  semilogxpanel(xlim=c(1e-6,1), xlab="Prey mass:predator mass",
                ylim=c(0,1), ylab="Preference")
  
  y = lapply(m, function(m) Phi(1/m, p, Delta=11))
  lines(m, y, lwd=thick, col=stdgrey)
  y = lapply(m, function(m) Phi(1/m, p, Delta=80))
  lines(m, y, lwd=thick)
  lines(m, phi(1/m,beta=p$beta,sigma=p$sigma), lty=dotted)
  
  legend("topleft", bty="n",
         legend=TeX(c("n=$\\infty, \\Delta$=1",
                  "n=10, $\\Delta$=11",
                  "n=6, $\\Delta$=80")),
         lty=c(dotted,solid,solid),
         lwd=c(1,thick,thick),
         col=c(black,stdgrey,black)
         )
  # sigma = function(Delta) {
  #   y = as.numeric( lapply(m, function(m) fTot(1,m,Delta)) )
  #   n0 = trapz(log(m),y)
  #   n1 = trapz(log(m), log(m)*y)/n0
  #   n2 = trapz(log(m), (log(m)-n1)^2*y)/n0
  #   return( sqrt(n2) )
  # }
  # 
  # semilogxpanel( xlim=Delta, xlab="Grid expansion factor $\\Delta$",
  #                ylim=c(1,2.5), ylab="Width $\\sigma$")
  # Delta = 10^seq(0,2,length.out=20)
  # y = lapply( Delta, sigma )
  # lines(Delta, y, lwd=thick)
  # lines(Delta, Delta/Delta*sigma, lty=dotted)
  # lines(c(10,10), c(1,3), lty=dotted)
}

plotGridtest = function() {
  n = c(6, 8, 10, 15, 25, 50)
  
  defaultplot()
  loglogpanel( xlim=c(5e-9,5000), xlab="Cell weight ($\\mu$gC)",
               ylim=c(2,200), ylab="Normalized biomass ($\\mu$gC/l)")
  
  textLegend = NULL
  for (i in 1:length(n)) {
    p = parametersChemostat( parameters(n[i]) )
    sim = simulate(p)
    Delta = p$m[2]/p$m[1]
    textLegend = c(textLegend, TeX(sprintf("n = %2.0f, $\\Delta$ = %2.1f", n[i],Delta)))
    
    Bnorm = sim$B/(sqrt(Delta)-1/sqrt(Delta))
    Bnorm[ Bnorm<1e-4 ] = NA
    lines( p$m, Bnorm )
    points( p$m, Bnorm, pch=i)
  }
  
  legend(x="topright", bty="n",
         legend=textLegend, #c("n=6","n=8","n=10","n=15","n=25","n=50"),
         pch = seq(1,6),
  )
}



