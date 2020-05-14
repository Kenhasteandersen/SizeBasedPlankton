source("modelWatercolumn.R")
source("basetools.R")

plotFunctionsWatercolumnDiffusion = function( p = parametersWatercolumn() ) {
  diff = 10^seq(-1,1, length.out=6)
  
  
  func = data.frame()
  for (i in 1:length(diff)) {
    p$diff = diff[i]
    
    sim = simulateWatercolumn(p)
    func = rbind( func, as.data.frame(calcFunctionsWatercolumn(sim) ))
  }
  func$diff = diff
  
  defaultplotvertical(nPanel=2)
  
  semilogxpanel(xlim=diff, ylim=c(0,0.05))
  lines(func$diff, func$Bpico, lwd=0.5)
  lines(func$diff, func$Bnano, lwd=1)
  lines(func$diff, func$Bmicro, lwd=1.5)
  lines(func$diff, func$Bpico+func$Bnano+func$Bmicro)
  
  semilogxpanel(xlim=diff, ylim=c(0,1000))
  lines(diff, func$prodCgross)
  lines(diff, func$prodCnet)
  lines(diff, func$prodHTL)
  
  return(func)
}
