library(deSolve)

phi = function(z, beta=500, sigma=1.3)
  exp( -(log(z/beta))^2/(2*sigma^2) )

calcrates = function(t,N,B) {
  JN = AN*m^(1/3)*N
  JL = AL*m^(2/3)*L
  for (i in 1:n) # loop over predators
    F[i] = sum(theta[i,]*B)
  JF = AF*m*F
  
  JCtot = JL+JF-Jresp
  JNtot = JN*rhoCN+JF # In units of carbon
  Jtot = pmin( JCtot, JNtot)
  
  f = (Jtot) / (Jtot + Jmax) # feeding level
  
  mortpred = rep(0,n)
  for (j in 1:n)  # loop over prey
    mortpred[j] = sum(theta[,j]*(1-f)*AF*B)
  
  dBdt = Jmax*f*B/m  - (mort + mortpred + mort2*B + mortHTL*(m>=mHTL))*B
  dNdt = d*(N0-N) - sum(f*JN*B/m)  # The uptake term is not quite right. It ignores Liebig
  
  return(list(m=m, N=N, B=B, dNdt=dNdt, dBdt=dBdt, JN=JN, JL=JL, JF=JF, Jtot=Jtot, f=f, mortpred=mortpred, mort=mort,
              totKilled = sum(JF*(1-f)*B/m), totEaten = sum(mortpred*B), totGrowth=sum(Jmax*f*B/m)))  
}

plotrates = function(t,N,B) {
  r = calcrates(t,N,B)
  m = r$m
  
  par(mfcol=c(2,1))

  plot(m, r$B, log="xy", type="b")
  lines(m,0*m,lty=3)
  
  plot(m, r$JN/m, log="x", type="l", col="blue", ylim=c(-1,1))
  lines(m, r$JL/m, col="green")
  lines(m, r$JF/m, col="red")
  lines(m, r$Jtot/m, lty=2)
  lines(m, Jmax/m*r$f)
  lines(m, r$dBdt/B,lwd=2)
  lines(m,0*m,lty=3)
  lines(m, -r$mortpred, col="red")
  lines(m, -r$mort, col="red", lty=2)
  lines(m, -Jresp/m, col="magenta")
  lines(m, -mort2*B, col="blue")
  lines(m, -mortHTL*(m>=mHTL), col="orange")
  
  
 # plot(m, r$dBdt/B, log="x", type="l", lwd=2, ylim=c(-0.2,1))
#  lines(m, r$mortpred, col="red")
#  lines(m, r$mort, col="red", lty=2)
#  lines(m,0*m,lty=3)
}

derivative = function(t,y,par) {
  N = y[1]
  B = y[2:(1+n)]
  
  rates = calcrates(t,N,B)
  
  return(list(c(rates$dNdt, rates$dBdt)))
}

n = 12; # No of groups
m = 10^seq(-7,0,length.out = n)  # mu g N
mHTL = m[n-2]

rhoCN = 106/16; # C:N ratio
#
# Clearance rates:
#
factor = 1.5^3*rhoCN*1e-6
AN = 2.5e-3*factor^(1/3)# l/day/mugN^1/3
AL = 72*factor^(2/3)#  mugC/day/(W^2 m^2)/mugN^(2/3)
AF = 1e-1*1.3e3*factor # l/d/mugN

par(mfcol=c(1,3))
plot(m, AN*m^(1/3), log="xy", type="b", xlab="mugN", ylab="l/day")
plot(m, AL*m^(2/3), log="xy", type="b", xlab="mugN", ylab="mugC/day/(W m^2)")
plot(m, AF*m, log="xy", type="b", xlab="mugN", ylab="l/day")
#
# Prey encounter
#
theta = matrix(nrow=n, ncol=n)
for (i in 1:n)
  theta[i,] = phi(m[i]/m)   # (predator, prey)


Jmax = 1 * m # mugC/day
Jresp = 0.1*Jmax
mort = 0*0.005*(Jmax/m) * m^(-1/4)
mort2 = 0.0015
mortHTL = 0.2

d = 1 # diffusion rate
N0 = 10 # Deep nutrient levels
B0 = rep(10,n)
L = 100

derivative(c(N0,B0),0)

nSave = 100
out = ode(c(N0,B0), seq(0, 1000, length.out = nSave), derivative)
ixB = 3:(n+2)

par(mfcol=c(2,1))
t = out[,"time"]
plot(t, out[,"1"], log="y", ylim=c(1e-5,1000), type="l", col="blue")
for (j in 1:n)
  lines(t, out[,j+2], lwd=j/3)

plot(m, out[nSave, 3:(n+2)], type="b", log="xy")

r=calcrates(t[nSave], out[nSave,2], out[nSave,ixB])
plotrates(t[nSave], out[nSave,2], out[nSave,ixB])
