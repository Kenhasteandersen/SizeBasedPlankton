require(latex2exp)
require(lattice)
source("basetools.R")
source("model.R")

# ===================================================
# Plots for documentation
# ===================================================

plotAll = function() {
  pdfplot("../AL.pdf", plotAL, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../AF.pdf", plotAF, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../AN.pdf", plotAN, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Mumax.pdf", plotMumax, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Rstar.pdf", plotRstar, width=doublewidth, height = height)
  
  plotSimulationExamples()
  
  pdf(file="../VaryLightAndDiffusion.pdf", width=doublewidth+0.5, height=height+.3)
  fontsize = trellis.par.get("fontsize")
  fontsize$text = 10
  trellis.par.set("fontsize",fontsize)
  plotVaryLightAndDiffusion()
  dev.off()
  
  pdfplot("../Functions.pdf", plotFunctions,n=25, width=doublewidth, height=2*height)
  
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
  points(data$C[!ixDiatom], data$A[!ixDiatom],pch=16, col="darkgreen")
  points(data$C[ixDiatom], data$A[ixDiatom],pch=16, col="darkgreen")
  points(data$C[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="darkgreen")
  
  C = 10^seq(log10(min(data$C)), log10(max(data$C)), length.out = 100)
  lines(C, exp(predict(fit, list(C=C))), lwd=2)
  #lines(C, exp(predict(fit_diatoms, list(C=C))), lwd=2, col="red")
  #lines(C, exp(predict(fit_no_diatoms, list(C=C))), lwd=3, col="blue")
  
  lines(C, 2*AL * C^(2/3), lty=dotted)
  lines(C, 0.026*C, lty=dotted)
  #lines(C, AL_diatoms * C^(2/3), col="red", lty=dotted)
  #lines(C, AL_no_diatoms * C^(2/3), col="blue", lty=dotted)
  
  # Camila's parameters:
  #cL=0.08733
  #AL=0.004864  
  #lines(C, cL*C * AL*C^(2/3) / ( cL*C + AL*C^(2/3) ) , col="orange")
  
  #legend(x="bottomright", bty="n",
  #       legend=c("Diatoms", "Other phototrophs", "Camila"),
  #       pch=c(16,16,NA),
  #       lwd=c(0,0,2),
  #       col=c("red","Blue","Orange"))
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
              ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/($\\mu$ mol photons m$^{-2}s^{-1})")
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
  loglogpanel(xlim = c(1e-9,1), ylim = A$mu_max,
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
  lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)), lwd=1)
  lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)) - parameters()$cR, lwd=2)
  lines(C, predict(fit, list(C=C)))
  
  # Camila
  #lines(C, 0.12*C^-0.25, col="orange", lwd=1)
  #lines(C, (0.12-0.03)*C^-0.25, col="orange", lwd=2)
  
  legend(x="bottomright", bty="n",
                  legend=c("Diatoms","Other phototrophs",
                           "Mixotrophs", "Heterotrophs",
                           "Bacteria","Model"),
                  lty=solid,
                  lwd=2,
                  col=c("darkgreen","green","blue","red","brown","black"))
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
  
  # Camila
  C = range(data$w[ixProtist])
  lines(C, 0.018*C, col="orange")
  
  
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
  for (i in taxons) {
    ix = dat$taxon == i
    points(C[ix], Anit[ix], pch=16, col=col)
    points(C[ix], Aamm[ix], pch=17, col=col)
    points(C[ix], Aphos[ix], pch=18, col=col)
    col = col + 1
  }
  #points(C, Aphos,pch=18)
  legend(x="bottomright",bty="n", 
         legend=taxons, pch=16, col=seq(1,length(taxons)))
  
  
  m = 10^seq(-7,1) # mu gC
  r = 1.5/2*(1e-6*m)^(1/3) # cm
  Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
  ANmax = 4*pi*Diff*r*1e-3 # L/day
  cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
  lines(m, ANmax, lwd=1, lty=dotted)  
  lines(m, parameters()$AN*m^(1/3), lwd=2)
  #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
  cat("alphaN = ", mean(parameters()$AN*m^(1/3)), "\n")
  
  V = 10^seq(-2,8,length.out = 100)
  alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
  lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
  
  # Camila's parameters:
  lines(C, 3.75e-5*C^(1/3), col="orange")
  
  #alphaN_Banas2011 = 2.6/0.1 / 14 / 6 # convert from uM N-> ugN, from gN->gC
  #lines(m , alphaN_Banas2011*m^(1-0.45), col="blue", lty=dashdotted)
  
  legend(x="topleft", bty="n",
         legend=c("Fit","Theoretical limit","Ward et al (2018)", "Camila"),
         lty=c(solid,dotted,dashed,solid),
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

plotRstar = function() {
  p = parameters()
  
  m = 10^seq(log10(0.3*min(p$m)), log10(max(p$m)), length.out = 1000)
  nu = pmin(1,p$c * m^(-1/3))
  mu = 0.4#c(0.2, 0.4, 0.6)
  
  tightaxes()
  defaultplot()
  loglogpanel(xlim=m, 
              ylim=c(0.005,1000),
              xlab = "Cell mass ($\\mu$gC)",
              ylab = "Limiting concentration ($\\mu$mol/l) or ($\\mu$E/m$^2$/s)")
  
  for (i in 1:length(mu)) {
    correct = mu[i]*(1-nu)/((1-nu)-mu[i])
    
    Nstar = m^0.667 /(p$AN*p$rhoCN) * correct
    Nstar[Nstar<0] = Inf
    lines(m, Nstar/14, lwd=2, col="blue")
    
    DOCstar = m^0.667/(p$AN) * correct
    DOCstar[DOCstar<0] = Inf
    lines(m, DOCstar/12, lwd=2, col="magenta")
    
    #Lstar = m^0.333 * p$alphaJ/(p$AL) * mu[i]/(p$alphaJ*(1-nu)-mu[i])
    #Lstar = p$alphaJ*(p$AL+p$cL*m^(1/3))*mu[i] / (p$AL*p$cL*(p$alphaJ-mu[i])*(1-nu))
    Lstar = m^(1/3)/(p$epsilonL*p$AL) *  1/(1-exp(-p$cL*m^(1/3))) * correct
    Lstar[Lstar<0] = Inf
    lines(m , Lstar, lwd=2, col="green")
  }
  
  legend(x="bottomright", bty="n",
         legend=c(TeX("$\\textit{N}^*$"), TeX("$$DOC^*$"), TeX("$$\\textit{L}^*$")),
         col=c("magenta","blue","green"),
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
  p = parameters()
  p$d = 0.0025
  p$L = 40
  p$mortHTL = 0.03
  
  sim = simulate(p)
  pdfplot("../spectrum1.pdf",FUN=plotSpectrum, sim=sim, width=doublewidth, height=1.5*height)
  pdfplot("../rates1.pdf", FUN=plotRates, sim=sim, width=doublewidth, height=1.5*height)
  #
  # Eutrophic
  #
  p = parameters()
  p$d = 0.08
  p$L = 20
  p$mortHTL = 0.1
  
  sim = simulate(p)
  pdfplot("../spectrum2.pdf",FUN=plotSpectrum, sim=sim, width=doublewidth, height=1.5*height)
  pdfplot("../rates2.pdf", FUN=plotRates, sim=sim, width=doublewidth, height=1.5*height)
}


plotVaryLightAndDiffusion = function() {
  library(lattice)
  library(latticeExtra)
  library(reshape)
  
  d = seq(0.002,.2,length.out=10) #10^seq(-2,log10(2),length.out = 10)
  L = seq(1,70,length.out=10)
  #
  # Simulations
  #
  df = data.frame(size=NULL, d=NULL, L=NULL, biomass=NULL)
  df2 = data.frame(func=NULL, d=NULL, L=NULL, z=NULL)
  #B = array(dim=c(3,length(d), length(L)))
  for (i in 1:length(d))
    for (j in 1:length(L)) {
      p = parameters()
      p$d = d[i]
      p$L = L[j]
      #p$mHTL = 0.008
      sim = simulate(p)
      func = calcFunctions(sim$p, sim$rates, sim$N, sim$B)
      
      df = rbind(df, data.frame(size="Pico",d=d[i],L=L[j],biomass=func$Bpico))
      df = rbind(df, data.frame(size="Nano",d=d[i],L=L[j],biomass=func$Bnano))
      df = rbind(df, data.frame(size="Micro",d=d[i],L=L[j],biomass=func$Bmicro))
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
                ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), xlab=TeX("Exchange rate (day$^{-1}$)"))
  
  plt2 = levelplot(z ~ d*L | func, data=df2,
                   aspect=1, region=TRUE,
                   col.regions = terrain.colors(100),
                   layout=c(3,1),
                   ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), 
                   xlab=TeX("Exchange rate (day$^{-1}$)"))
  
  return(plt)
}

plotFunctions = function(L=c(20,50), n=10) {
  
  panelsFunctions = function(L=c(18,40), n=10, yaxis=TRUE) {
    d = 10^seq(-3,log10(.2),length.out = n) #seq(0.02,2,length.out=n) #
    
    p = parameters()
    p$L = L
    
    F = data.frame()
    for (i in 1:length(d)) {
      p$d = d[i]
      sim = simulate(p)
      func = calcFunctions(sim$p, sim$rates, sim$N, sim$B)
      func$d = d[i]
      func$N = sim$N
      func$DOC = sim$DOC
      
      F = rbind(F, as.data.frame(func))
    }
    
    # Biomass
    B = F$Bpico+F$Bnano+F$Bmicro
    defaultpanel(xlim=d, ylim=c(0,0.15), xaxis=FALSE, yaxis=yaxis, bty="l",
                 ylab="Biomass ($\\mu$gC/l)")
    lines(F$d, B, lwd=2)
    
    # N
    semilogypanel(xlim=d, ylim=c(0.0001,10), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="N ($\\mu$mol N/l)")
    lines(F$d, F$N/14, lwd=2)
    
    #DOC
    semilogypanel(xlim=d, ylim=c(1,500), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="DOC  ($\\mu$mol C/l)")
    lines(F$d, 1000*F$DOC/12, lwd=2)
    
    #Production:
    semilogypanel(xlim=d, ylim=c(1,2000), xaxis=FALSE, yaxis=yaxis, bty="l",
                  ylab="Prod. (gC/m$^2$/yr)")
    lines(F$d, F$prodCgross)
    lines(F$d, F$prodCnet, col="blue")
    lines(F$d, F$prodNew, col="green")
    lines(F$d, F$prodHTL, col="red")
    
    # Eff
    defaultpanel(xlim=d, ylim=c(0,1), bty="l",
                 ylab="$\\epsilon_{HTL}$",
                 xlab="Exchange rate (day$^{-1}$)",
                 yaxis=yaxis)
    lines(F$d, F$effHTL, lwd=2)
  }
  
  defaultplotvertical(mfcol=c(5,2))
  
  panelsFunctions(L=L[1],n=20)
  panelsFunctions(L=L[2], yaxis=FALSE,n=20)
}
