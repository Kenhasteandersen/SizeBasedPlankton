// Build: g++ -O3 -shared model.cpp -o model.so  

#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h> 

typedef std::vector<double> type_mass;  // type for bins

struct plankton {
  int n; // No of size classes
  type_mass m;
  double rhoCN;
  double epsilonL;
  double epsilonF;
  type_mass ANm, ALm, AFm;
  type_mass Jmax, Jresp, JFmax;
  type_mass Jloss_passive;
  
  std::vector< std::vector<double> > theta;
  
  type_mass mort;
  double mort2;
  type_mass mHTL;
  
  double remin, remin2;
  double cLeakage;
};

struct typeRates {
  type_mass JN, JDOC, JL, JF, F, JFreal;
  type_mass JNtot, JLreal, JCtot, Jtot;
  type_mass JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig;
  type_mass JNloss, JCloss;
  type_mass mortpred;
};

plankton p;
typeRates rates;
type_mass B, ANmT, JmaxT, JFmaxT, JrespT;
const int idxN = 0;
const int idxDOC = 1;
const int idxB = 2;

inline double min(double a, double b) {
  return (a < b) ? a : b;
};

inline double max(double a, double b) {
  return (a < b) ? b : a;
};


void printRate(type_mass J) {
  for (int i=0; i<p.n; i++) 
    std::cout << J[i]/p.m[i] << ", ";
  std::cout << "\n";
}

void printRates() {
  std::cout << "jN: ";
  printRate(rates.JN);
  
  std::cout << "jF: ";
  printRate(rates.JF);

  std::cout << "jFreal: ";
  printRate(rates.JFreal);
  
  std::cout << "jL: ";
  printRate(rates.JLreal);
  
  std::cout << "jLreal: ";
  printRate(rates.JLreal);
  
  std::cout << "jCtot: ";
  printRate(rates.JCtot);
  
  std::cout << "jNtot: ";
  printRate(rates.JNtot);
  
  std::cout << "jTot: ";
  printRate(rates.Jtot);
}

/* ===============================================================
 * Stuff for cell model:
 */
extern "C" void setParameters(
    int& _n, 
    double* _m,
    double& _rhoCN, 
    double& _epsilonL, 
    double& _epsilonF,
    double* _ANm, 
    double* _ALm, 
    double* _AFm,
    double* _Jmax, 
    double* _JFmax,
    double* _Jresp,
    double* _Jloss_passive,
    double* _theta,
    double* _mort,
    double& _mort2,
    double* _mortHTL,
    double& _remin,
    double& _remin2,
    double& _cLeakage
  ) {
  
  p.n = _n;
  p.rhoCN = _rhoCN;
  p.epsilonL = _epsilonL;
  p.epsilonF = _epsilonF;
  p.mort2 = _mort2;
  p.remin = _remin;
  p.remin2 = _remin2;
  p.cLeakage = _cLeakage;
  
  p.m.resize(p.n);
  p.ANm.resize(p.n);
  p.ALm.resize(p.n);
  p.AFm.resize(p.n);
  p.Jmax.resize(p.n);
  p.JFmax.resize(p.n);
  p.Jresp.resize(p.n);
  p.Jloss_passive.resize(p.n);
  p.mort.resize(p.n);
  p.mHTL.resize(p.n);
  
  p.theta.resize(p.n);
  for (int i=0; i<p.n; i++) {
    p.m[i] = _m[i];
    p.ANm[i] = _ANm[i];
    p.ALm[i] = _ALm[i];
    p.AFm[i] = _AFm[i];
    p.Jmax[i] = _Jmax[i];
    p.JFmax[i] = _JFmax[i];
    p.Jresp[i] = _Jresp[i];
    p.Jloss_passive[i] = _Jloss_passive[i];
    p.mort[i] = _mort[i];
    p.mHTL[i] = _mortHTL[i];
    
    p.theta[i].resize(p.n);
    for (int j=0; j<p.n; j++)
      p.theta[i][j] = _theta[i*p.n+j];
  };
  //
  // Init private temps:
  //
  B.resize(p.n);
  ANmT.resize(p.n);
  JmaxT.resize(p.n);
  JFmaxT.resize(p.n);
  JrespT.resize(p.n);
  //
  //  Init rates:
  //
  rates.JN.resize(p.n);
  rates.JDOC.resize(p.n);
  rates.JL.resize(p.n);
  rates.JF.resize(p.n);
  rates.JFreal.resize(p.n);
  rates.F.resize(p.n);
  rates.JNtot.resize(p.n);
  rates.JLreal.resize(p.n);
  rates.JCtot.resize(p.n);
  rates.Jtot.resize(p.n);
  
  rates.JCloss_feeding.resize(p.n);
  rates.JCloss_photouptake.resize(p.n);
  rates.JNlossLiebig.resize(p.n);
  rates.JClossLiebig.resize(p.n);
  
  rates.JNloss.resize(p.n);
  rates.JCloss.resize(p.n);
  
  rates.mortpred.resize(p.n);
};

inline double fTemp(double Q10, double T) {
  return pow(Q10, T/10-1);
};

void calcRates(const double& T, const double& L, const double* u, 
               double* dudt, const double gammaN, const double gammaDOC) {
  int i;
  double f;
  
  for (i=0; i<p.n; i++)
    B[i] = (0 < u[idxB+i]) ? u[idxB+i] : 0;
  //
  // Temperature corrections:
  //
  double f15 = fTemp(1.5, T);
  double f20 = fTemp(2.0, T);
  for (i=0; i<p.n; i++) {
    ANmT[i] = p.ANm[i]*f15;
    JmaxT[i] = p.Jmax[i]*f20;
    JrespT[i] = p.Jresp[i]*f20;
    JFmaxT[i] = p.JFmax[i]*f20;
  }
  //
  // Uptakes
  //
  for (i=0; i<p.n; i++) {
    /*
     * Uptakes:
     */ 
    
    // Nutrients:
    if (u[idxN]<=0)
      rates.JN[i] = 0;
    else
      rates.JN[i] = gammaN * ANmT[i]*u[idxN]*p.rhoCN;
    
    // DOC:
    if (u[idxDOC]<=0)
      rates.JDOC[i] = 0;
    else
      rates.JDOC[i] = gammaDOC * ANmT[i]*u[idxDOC];
    
    // Light:
    rates.JL[i] = p.epsilonL * p.ALm[i]*L;
    
    // Feeding:
    rates.F[i] = 0;
    for (int j = 0; j<p.n; j++) {
      rates.F[i] += p.theta[j][i]*B[j];
      //std::cout << rates.F[i] << ";";
    }
    rates.JF[i] = p.epsilonF * JFmaxT[i] * p.AFm[i]*rates.F[i] 
      / (p.AFm[i]*rates.F[i] + JFmaxT[i]);
    
    // Combined N and C fluxes:
    rates.JNtot[i] = rates.JN[i] + rates.JF[i] - p.Jloss_passive[i];
    rates.JCtot[i] = rates.JL[i] + rates.JF[i] + rates.JDOC[i]  
      - JrespT[i] - p.Jloss_passive[i];

    // Synthesis:
    rates.Jtot[i] = min( rates.JNtot[i], rates.JCtot[i] ); 
    f = rates.Jtot[i]/(rates.Jtot[i] + JmaxT[i]);
      
    //If synthesis-limited then down-regulate feeding:
    rates.JFreal[i] = max(0, rates.JF[i] - (rates.Jtot[i]-f*JmaxT[i]));
    rates.Jtot[i] = f*JmaxT[i];
    rates.JLreal[i] = rates.JL[i] 
    - min((rates.JCtot[i] - (rates.JF[i]-rates.JFreal[i])-f*JmaxT[i]), rates.JL[i]);
        
    // Actual uptakes:
    rates.JCtot[i] = rates.JLreal[i] + rates.JDOC[i] 
    + rates.JFreal[i] - p.Jresp[i] - p.Jloss_passive[i];
    rates.JNtot[i] = rates.JN[i] + rates.JFreal[i] - p.Jloss_passive[i];
  
    // Losses:
    rates.JCloss_feeding[i] = (1-p.epsilonF)/p.epsilonF*rates.JFreal[i];
    rates.JCloss_photouptake[i] = (1-p.epsilonL)/p.epsilonL*rates.JLreal[i];
    rates.JNlossLiebig[i] = max(0, rates.JNtot[i]-rates.Jtot[i]);
    rates.JClossLiebig[i] = max(0, rates.JCtot[i]-rates.Jtot[i]); 

    rates.JNloss[i] = rates.JCloss_feeding[i]
      +rates.JNlossLiebig[i]
      +p.Jloss_passive[i];
    rates.JCloss[i] = rates.JCloss_feeding[i] 
      + rates.JCloss_photouptake[i] 
      + rates.JClossLiebig[i]
      + p.Jloss_passive[i];
  }
  //
  // Mortality:
  //
  for (i=0; i<p.n; i++) {
    rates.mortpred[i] = 0;
    for (int j = 0; j<p.n; j++) {
      rates.mortpred[i] += p.theta[i][j] * rates.JFreal[j]*B[j]/(p.epsilonF*p.m[j]*rates.F[j]);
    }
  }
  //
  // Reaction terms:
  //
  dudt[idxN] = 0;
  dudt[idxDOC] = 0;
  for (i=0; i<p.n; i++) {
    double mortloss;
    mortloss = B[i]*((1-p.remin2)*p.mort2*B[i] + p.mHTL[i]);
    dudt[idxN] += 
      ((-rates.JN[i]
      + rates.JNloss[i])*B[i]/p.m[i] 
      + p.remin2*p.mort2*B[i]*B[i]
      + p.remin*mortloss)/p.rhoCN;
      
    dudt[idxDOC] += 
      (-rates.JDOC[i] 
      + rates.JCloss[i])*B[i]/p.m[i] 
      + p.remin2*p.mort2*B[i]*B[i]
      + p.remin*mortloss;
      
    dudt[idxB+i] = (rates.Jtot[i]/p.m[i]  
      - (p.mort[i] 
        + rates.mortpred[i] 
        + p.mort2*B[i] 
        + p.mHTL[i]))*B[i];

  }
  //printRates();
};

extern "C" void calcRates(const double& T, const double& L, const double* u, double* dudt) {
  calcRates(T,L,u,dudt,1,1);
}

void calcRates(const double& T, const double& L, const double* u, double* dudt, const double dt) {
  calcRates(T, L, u, dudt);
  
  double gammaN = 1;
  double gammaDOC = 1;

  if (u[idxN]+dudt[idxN]*dt<0)
    gammaN = min(gammaN, -u[idxN]/(dudt[idxN]*dt));
  
  if ((u[idxDOC]+dudt[idxDOC]*dt)<0) {
    //double JDOCtot = 0;
    //for (int j = 0; j<p.n; j++)
    //  JDOCtot += rates.JDOC[j];
    min(gammaDOC, -u[idxDOC]/(dudt[idxDOC]*dt));
  }
    
  if (gammaN<1 | gammaDOC<1) {
    calcRates(T,L,u,dudt,gammaN,gammaDOC);
    std::cout << gammaN << "," << gammaDOC << "\n";
  }
};

/* ===============================================================
 * Stuff for Chemostat model:
 */
extern "C" void derivativeChemostat(const double& L, const double& T, const double& d, const double& N0,
                                   const double* u, double* dudt) {
  calcRates(T, L, u, dudt);
  
  dudt[idxN]     += d*(N0-u[idxN]);
  dudt[idxDOC]   += d*(0-u[idxDOC]);
  
  for (int i=0; i<p.n; i++)
    dudt[idxB+i] += d*(0-u[idxB+i]);
};
/* ===============================================================
 * Stuff for watercolumn model:
 */

inline double calcLight(double L0, double x) {
  double k = 0.1;
  return L0*exp(-k*x);
};

void solveTridiag(std::vector<double>&  a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& s, 
                  const double* dudt, const double dt, const int step, const double* uprev,
                  std::vector<double>& cprime, std::vector<double>& sprime, double* u) {
  int i;
  int n = cprime.size();
  
  cprime[0] = c[0]/b[0];
  sprime[0] = (s[0]+dudt[0]*dt+uprev[0])/b[0];
  for (i=1; i<n; i++) {
    cprime[i] = c[i]/(b[i]-a[i]*cprime[i-1]);
    sprime[i] = (s[i]+dudt[i*step]*dt+uprev[i*step]-a[i]*sprime[i-1]) / (b[i]-a[i]*cprime[i-1]);
  }
  u[(n-1)*step] = sprime[(n-1)*step];
  for (i=n-1; i>=0; i--)
    u[i*step] = sprime[i] - cprime[i]*u[(i+1)*step];
};

void printu(const double* u, const int nGrid) {
  std::cout << "u:\n";
  for (int j=0; j<nGrid; j++) {
    std::cout << j << ": ";
    for (int k=0; k<(p.n+2); k++)
      std::cout << u[j*(p.n+2)+k] << ",";
    std::cout << "\n";
  }
  std::cout << "------------\n";
}

extern "C" void simulateWaterColumnFixed(const double& L0, const double& T, 
                                        const double& Diff, const double& N0,
                                        const double& tEnd, const double& dt,
                                        const int& nGrid, const double* x, double* u) {
  /*
   * Set up matrices
   */
  int i,j,k;
  double dx = x[1]-x[0];
  double Dif = dt/(dx*dx)*Diff;
  
  std::vector<double> aN, bN, cN, sN;
  std::vector<double> aDOC, bDOC, cDOC, sDOC;
  std::vector<double> aPhyto, bPhyto, cPhyto, sPhyto;
  std::vector<double> cprime, sprime;
  
  aN.resize(nGrid);
  bN.resize(nGrid);
  cN.resize(nGrid);
  sN.resize(nGrid);
  aDOC.resize(nGrid);
  bDOC.resize(nGrid);
  cDOC.resize(nGrid);
  sDOC.resize(nGrid);
  aPhyto.resize(nGrid);
  bPhyto.resize(nGrid);
  cPhyto.resize(nGrid);
  sPhyto.resize(nGrid);
  
  cprime.resize(nGrid);
  sprime.resize(nGrid);
  
  for (i=0; i<nGrid; i++) {
    aN[i] = -Dif;
    bN[i] = 1 + 2*Dif;
    cN[i] = -Dif;
    sN[i] = 0;
  }
  bN[0] = 1 + Dif;
  sN[nGrid-1] = -cN[nGrid-1]*N0;
  
  aDOC = aN;
  bDOC = bN;
  cDOC = cN;
  sDOC = sN;
  sDOC[nGrid-1] = 0;
  
  double adv = 0*dt/dx;
  aPhyto = aN;
  bPhyto = bN; // + adv;
  cPhyto = cN;
  sPhyto = sN;
  sPhyto[nGrid-1] = 0;
  /*
   * Initialize
   */
  double *dudt;
  dudt = (double *) malloc((p.n+2)*nGrid*sizeof(double *));
  
  /*
   * Iterate
   */
  for (i=1; i<tEnd/dt; i++) {
    // Calc reaction:
    for (j=0; j<nGrid; j++) {
      calcRates(T, calcLight(L0, x[j]), &u[j*(p.n+2) + (i-1)*(p.n+2)*nGrid], &dudt[j*(p.n+2)], dt);
    }
    // Invert:
    solveTridiag(aN, bN, cN, sN, 
                 &dudt[idxN], dt, p.n+2,
                 &u[nGrid*(p.n+2)*(i-1)+idxN],
                 cprime, sprime, &u[nGrid*(p.n+2)*i+idxN]);
    solveTridiag(aDOC, bDOC, cDOC, sDOC, 
                 &dudt[idxDOC], dt, p.n+2,
                 &u[nGrid*(p.n+2)*(i-1)+idxDOC],
                 cprime, sprime, &u[nGrid*(p.n+2)*i+idxDOC]);
    for (k=0; k<p.n; k++)
      solveTridiag(aPhyto, bPhyto, cPhyto, sPhyto, 
                   &dudt[idxB+k], dt, p.n+2,
                   &u[nGrid*(p.n+2)*(i-1)+idxB+k],
                   cprime, sprime, &u[nGrid*(p.n+2)*i+idxB+k]);
  }  
  
}

int main()
{
  return 0;
};
