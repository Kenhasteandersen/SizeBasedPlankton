source("model.R")
source("basetools.R")
library(tictoc)

parametersWatercolumn = function(p = parameters()) {
  p$tEnd = 2*365
  p$dt = 0.01
  p$depth = 60
  p$Lsurface = 100
  p$nGrid = 50
  p$k = 0.1 # Damping of light by water
  p$diff = 1
  
  return(p)
}

simulateWatercolumn = function(p=parametersWatercolumn(), 
                               tEnd=p$tEnd, dt=p$dt, nGrid=p$nGrid, nTout=110) {
  # Load C engine
  if (!is.loaded(paste(fileLibrary,".so",sep="")))
    dyn.load(paste(fileLibrary,".so",sep=""))
  # Constants:
  idxN = 1
  idxDOC = 2
  idxB = 3:(2+p$n)
  #
  # Set parameters in C engine:
  #
  mortHTL = 0*p$m + p$mortHTL*(p$m>p$mHTL)
  dummy = .C("setParameters", as.integer(p$n), 
             p$m, p$rhoCN, p$epsilonL, p$epsilonF,
             p$ANm, p$ALm, p$AFm, p$Jmax, p$Jresp, p$theta,
             p$mort, p$mort2, mortHTL, p$remin);
  #
  # Set up grid and solution matrix:
  #
  x = as.double(seq(0, p$depth, length.out = nGrid) )
  
  u = array(0, dim=c(2+p$n, nGrid, tEnd/dt+2));
  # Initial conditions:
  u[idxN,,1] = p$N0;
  u[idxDOC,,1] = p$DOC0;
  u[idxB,,1] = p$B0;
  u[idxB, x>50, 1] = 1e-8;
  #
  # Simulate
  tictoc::tic()
  sim = .C("simulateWaterColumnFixed", 
           L0=as.double(p$L),
           T=as.double(p$T),
           Diff=as.double(p$diff),
           N0 = as.double(p$N0),
           tEnd=as.double(tEnd),
           dt=as.double(dt),
           nGrid=as.integer(nGrid),
           x=x,
           u=u)
  tictoc::toc()
  #
  # Extract results:
  #
  idxT = seq(1,tEnd/dt-2,length.out = nTout)
  idxX = seq(p$nGrid, 1, by=-1)
  
  sim$N = t(sim$u[idxN, idxX, idxT]);
  sim$DOC = t(sim$u[idxDOC, idxX, idxT])
  sim$B = sim$u[idxB, idxX, idxT]
  d = calcESD(p$m)
  sim$Bpico = t(colSums(sim$B[d<2,,],1))
  sim$Bnano = t(colSums(sim$B[d>=2 & d<20,,],1))
  sim$Bmicro = t(colSums(sim$B[d>=20,,],1))
  
  sim$t = seq(0, tEnd, by=dt)[idxT]
  sim$x = -x[seq(nGrid,1,by=-1)]
  
  sim$p = p
  return(sim)
}

plotWatercolumnTime = function(sim) {
  library(fields)
  #
  # Plot results
  #
  defaultplotvertical(nPanels = 5)
  image(sim$t, sim$x, sim$N, ylab="N", xlab="time", col=topo.colors(20))
  mtext(side=left, line=1, TeX("N"), cex=par()$cex)
  image(sim$t, sim$x, sim$DOC, col=topo.colors(20), 
        zlim=c(0,max(sim$DOC[(length(s$t/2)):length(s$t),])))
  mtext(side=left, line=1, TeX("DOC"), cex=par()$cex)
  
  zlim = c(0, max( (sim$Bpico+sim$Bnano+sim$Bmicro)[(length(s$t/2)):length(s$t),]))
  image(sim$t, sim$x, sim$Bpico, ylab="Bpico", col=topo.colors(20), zlim=zlim)
  mtext(side=left, line=1, TeX("Bpico"), cex=par()$cex)
  image(sim$t, sim$x, sim$Bnano, ylab="Bnano", col=topo.colors(20), zlim=zlim)
  mtext(side=left, line=1, TeX("Bnano"), cex=par()$cex)
  image(sim$t, sim$x, sim$Bmicro, xlab="time (days)", 
        ylab="Bmicro", col=topo.colors(20), zlim=zlim)
  mtext(side=left, line=1, TeX("Bmicro"), cex=par()$cex)
  mtext(side=bottom, line=1, TeX("Time (days)"), cex=par()$cex)
  
}

plotWatercolumn = function(sim, idx = length(sim$t)) {
  p = sim$p
  defaultplot(mfcol=c(1,2))
  #
  # Rates:
  #
  defaultpanel(xlim=c(0,80), xlab="Concentration ($\\mu$g C/l)",
               ylim=range(sim$x), ylab="Depth (m)")
  tightaxes()  
  lines(sim$N[idx,]/10, sim$x, col="blue", lwd=thick)
  lines(sim$DOC[idx,]*100, sim$x, col="magenta", lwd=thick)
  lines(sim$Bpico[idx,], sim$x, col="black", lwd=1)
  lines(sim$Bnano[idx,], sim$x, col="black", lwd=2)
  lines(sim$Bmicro[idx,], sim$x, col="black", lwd=3)
  
  title("\nConcentrations")
  legend(x="bottom", bty="n",
         legend=c("Nutrients/10", "DOC*100", "Pico", "Nano", "Micro"),
         lwd=c(thick, thick, 1,2,3),
         col=c("blue","magenta","black","black","black"))
  #
  # Carbon strategy:
  #
  rates=list()
  col = matrix(0,length(sim$x),p$n)
  pp = p
  for (j in 1:length(sim$x)) {
    rates[[j]] = calcRates(sim$t[idx], p$Lsurface*exp(p$k*sim$x[j]), sim$N[idx,j], sim$DOC[idx,j],
                           sim$B[,j,idx],pp)
    for (k in 1:length(p$m))
      col[length(sim$x)-j+1,k] = rgb(min(1,max(0,3*(rates[[j]]$JF/p$m)[k])), 
                                     min(1, max(0,3*(rates[[j]]$JLreal/p$m)[k])),
                                     min(1, max(0,3*(rates[[j]]$JDOC/p$m)[k])))
  }
  
  m = log10(p$m)
  dm = 0.5*diff(m)[1]
  defaultpanel(xlim=m, xlab="Cell mass (log10($\\mu$gC))",
               ylim=sim$x, ylab="")
  tightaxes()
  rasterImage(as.raster(col), min(m)-dm, min(sim$x), max(m)+dm, max(sim$x),
              interpolate = TRUE)
  
  title("\nCarbon strategy", col.main="white")
  legend(x="bottom", bty="n", cex=cex,
         legend=c("DOC", "Photoharvesting","Phagotrophy"),
         text.col="white",
         fill=c("blue", "green","red","transparent"),
         border=c("white","white","white","transparent"),
         lwd = c(0,0,0,3),
         col=c(NA,NA,NA,1))
}

calcFunctionsWatercolumn = function(sim) {
  param = sim$p
  idx = length(sim$t)
  conversion  = 365*1e-6*1000 # Convert to gC/yr/m2
  dx = diff(sim$x)
  dx = c(dx, dx[length(dx)])
  
  prodCgross = 0
  prodCnet = 0
  prodHTL = 0
  resp = 0
  Bpico = 0
  Bnano = 0
  Bmicro = 0

  p = param
  rates = NULL
  for (i in 1:param$nGrid) {
    # Calc rates throughout the water column
    L = param$Lsurface*exp(-p$k*sim$x[i])
    rates[[i]] = calcRates(0, L, sim$N[idx,i], sim$DOC[idx,i], sim$B[,i,idx],p)
    #
    # Primary production (carbon fixed)
    #
    prodCgross = prodCgross + 
      conversion * sum(rates[[i]]$JLreal*sim$B[,i,idx]/param$m)/param$epsilonL*dx[i]
    prodCnet = prodCnet + 
      conversion * sum( (rates[[i]]$JLreal)*sim$B[,i,idx]/param$m )*dx[i]
    #
    # Loss to HTL:
    #
    prodHTL = prodHTL + 
      conversion * sum( param$mortHTL*(param$m>=param$mHTL)*sim$B[,i,idx] )*dx[i]
    #
    # Loss to depth:
    #
    #prodSeq = 0#conversion * r$POMgeneration * (1-epsilonPOM)
    #
    # Respiration
    #
    resp = resp + conversion * sum( param$Jresp*sim$B[,i,idx]/param$m )*dx[i]
    #
    # Biomasses:
    #
    conversion2 = 1/10*100*1e-6  # Convert to gC/m2
    d = calcESD(param$m)
    Bpico = Bpico + conversion2 * sum( sim$B[d < 2,i,idx] )*dx[i]
    Bnano = Bnano + conversion2 * sum( sim$B[d>=2 & d <20,i,idx] )*dx[i]
    Bmicro = Bmicro + conversion2 * sum( sim$B[d>=20,i,idx])*dx[i]
  }
  #
  # Efficiencies:
  #
  effHTL = prodHTL/prodCgross*param$epsilonL # NOTE DIFFERENT DEFINITION THAN IN CHEMOSTAT
  
  return(list(
    prodCgross = prodCgross,
    prodCnet = prodCnet,
    prodHTL = prodHTL,
    resp = resp,
    effHTL = effHTL,
    Bpico=Bpico, Bnano=Bnano, Bmicro=Bmicro
  ))
}

#
# Plot functions:
#
plotFunctionsWatercolumn <- function(sim) {
  # Get the func value from the previous call:
  oldfunc = attr(plotFunctionsWatercolumn, "oldfunc")
  if (is.null(oldfunc))
    oldfunc = c(0,0,0)
  
  func = calcFunctionsWatercolumn(sim)
  fnc = c(func$prodCgross, func$prodCnet, func$prodHTL)
  attr(plotFunctionsWatercolumn, "oldfunc") <<- fnc
  heights = matrix(c(fnc, oldfunc), nrow=2, byrow = TRUE)
  
  par(mar=c(5,12,4,2))
  
  barplot(height=heights,
          names.arg = c("Gross PP", "Net PP", "HTL"),
          xlab = TeX("Production (gC/m$^2$/yr)"),
          beside=TRUE, col=c("black","grey"),
          horiz=TRUE, las=1,
          border=NA)
  legend("topright",
         c("This simulation","Previous simulation"),
         fill=c("black","grey"),
         bty="n")
}

baserunWatercolumn = function(p=parametersWatercolumn()) {
  sim = simulateWatercolumn(p)
  plotWatercolumn(sim)
  
  return(sim)
}

#
# Find cpp-code:
#
if (file.exists("../Cpp/model.cpp")) {
  fileLibrary = "../Cpp/model"
} else {
  if (file.exists("model.cpp")) {
    fileLibrary = "model"
  } else {
    stop("Did not find model cpp source code.")
}}
#
# Compile:
#
if (!file.exists(paste(fileLibrary,'.so',sep=""))) {
  if (system2("g++", 
           paste("-O3 -fPIC -shared ",fileLibrary,".cpp -o ",
               fileLibrary,".so",sep=""))==1)  {
    stop("Cannot compile cpp engine")
  }
}



