source("model.R")
source("basetools.R")

simulateWatercolumn = function(p=parameters(), 
                               tEnd=p$tEnd, dt=p$dt, nGrid=p$nGrid, nTout=110) {
  # Load C engine
  if (!is.loaded("../Cpp/model.so"))
    dyn.load("../Cpp/model.so")
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
             Diff=as.double(p$d),
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
  return(sim)
}

plotWatercolumnTime = function(p, sim) {
  library(fields)
  #
  # Plot results
  #
  defaultplotvertical(nPanels = 5)
  image.plot(sim$t, sim$x, sim$N)
  image.plot(sim$t, sim$x, sim$DOC)
  image.plot(sim$t, sim$x, sim$Bpico)
  image.plot(sim$t, sim$x, sim$Bnano)
  image.plot(sim$t, sim$x, sim$Bmicro, xlab="time (days)")
}

plotWatercolumn = function(p, sim, idx = length(sim$t)) {
  defaultplot(mfcol=c(1,2))
  
  defaultpanel(xlim=c(0,150), xlab="mu g C/?",
               ylim=range(sim$x), ylab="depth(m)")
  tightaxes()  
  lines(sim$N[idx,], sim$x, col="blue", lwd=thick)
  lines(sim$DOC[idx,], sim$x, col="magenta", lwd=thick)
  lines(sim$Bpico[idx,], sim$x, col="black", lwd=1)
  lines(sim$Bnano[idx,], sim$x, col="black", lwd=2)
  lines(sim$Bmicro[idx,], sim$x, col="black", lwd=3)

  rates=list()
  col = matrix(0,length(sim$x),p$n)
  pp = p
  for (j in 1:length(sim$x)) {
    pp$L = p$L*exp(0.1*sim$x[j])
    rates[[j]] = calcRates(sim$t[idx], sim$N[idx,j], sim$DOC[idx,j],
                         sim$B[,j,idx],pp)
    for (k in 1:length(p$m))
      col[length(sim$x)-j+1,k] = rgb(min(1,max(0,3*(rates[[j]]$JF/p$m)[k])), 
                  min(1, max(0,3*(rates[[j]]$JLreal/p$m)[k])),
                  min(1, max(0,3*(rates[[j]]$JDOC/p$m)[k])))
  }
  
  m = log10(p$m)
  dm = 0.5*diff(m)[1]
  defaultpanel(xlim=m, xlab="size, log10(mu gC)",
               ylim=sim$x, ylab="")
  tightaxes()
  rasterImage(as.raster(col), min(m)-dm, min(sim$x), max(m)+dm, max(sim$x),
              interpolate = TRUE)
}

baserunWatercolumn = function(p=parameters()) {
  p$L=300
  p$tEnd = 200
  p$dt = 0.1
  p$nGrid = 100
  p$depth = 100
  
  sim = simulateWatercolumn(p)
  plotWatercolumn(p, sim)
  
  return(sim)
}
