source("modelWatercolumn.R")
source("basetools.R")

plotAllWatercolumn = function() {
  pdfplot("../FunctionsWatercolumnDiffusion.pdf", plotFunctionsWatercolumnDiffusion, width = singlewidth, height=2*height)
}

plotFunctionsWatercolumnDiffusion = function( p = parametersWatercolumn(), n=10 ) {

  
  diff = 10^seq(-1,log10(1000), length.out=n)
  
  func = data.frame()
  for (i in 1:length(diff)) {
    p$diff = diff[i]
    
    sim = simulateWatercolumn(p)
    plotWatercolumn(sim)
    func = rbind( func, as.data.frame(calcFunctionsWatercolumn(sim) ))
  }
  func$diff = diff
  
  defaultplotvertical(nPanel=3)
  
  loglogpanel(xlim=diff, xaxis=FALSE,
                ylim=c(0.001,0.2), ylab="gC/m$^2$")
  lines(func$diff, func$Bpico, lwd=0.5, col=stdgrey)
  lines(func$diff, func$Bnano, lwd=1, col=stdgrey)
  lines(func$diff, func$Bmicro, lwd=1.5, col=stdgrey)
  lines(func$diff, func$Bpico+func$Bnano+func$Bmicro, lwd=thick)
  lines(diff, 100*func$BDOC, lwd=thick, col="magenta")
  
  loglogpanel(xlim=diff, xaxis=FALSE,
                ylim=c(1,10000), ylab="gC/m$^2/yr")
  lines(diff, func$prodCgross)
  lines(diff, func$prodCnet, col="blue")
  lines(diff, func$prodHTL,col="red")
  
  semilogxpanel(xlim=diff, xlab="Diffusion (m$^2$/day)",
                ylim=c(0,0.25), ylab="$\\epsilon_{HTL}$")
  lines(diff, func$effHTL, lwd=thick)
  
  return(func)
}
