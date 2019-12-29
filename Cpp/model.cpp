// Build: g++ -O3 -shared model.cpp -o model.so  
  
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>

typedef std::vector<double> type_mass;  // type for bins

struct plankton {
  int n; // No of size classes
  type_mass m;
  double rhoCN;
  double epsilonL;
  double epsilonF;
  type_mass ANm, ALm, AFm;
  type_mass Jmax, Jresp;

  std::vector< std::vector<double> > theta;

  type_mass mort;
  double mort2;
  type_mass mHTL;

  double remin;
};

struct typeRates {
  type_mass JN, JDOC, JL, JF, F;
  type_mass JNtot, JLreal, JCtot, Jtot;
  type_mass JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig;
  type_mass JNloss, JCloss;
  type_mass mortpred;
};

plankton p;
typeRates rates;
type_mass B, ANmT, JmaxT, JrespT;
const int idxN = 0;
const int idxDOC = 1;
const int idxB = 2;

inline double min(double a, double b) {
  return (a < b) ? a : b;
};

inline double max(double a, double b) {
  return (a < b) ? b : a;
};

extern "C" void setParameters(int& _n, double* _m,
  double& _rhoCN, double& _epsilonL, double& _epsilonF,
  double* _ANm, double* _ALm, double* _AFm,
  double* _Jmax, double* _Jresp,
  double* _theta,
  double* _mort,
  double& _mort2,
  double* _mortHTL,
  double& _remin) {

  p.n = _n;
  p.rhoCN = _rhoCN;
  p.epsilonL = _epsilonL;
  p.epsilonF = _epsilonF;
  p.mort2 = _mort2;
  p.remin = _remin;

  p.m.resize(p.n);
  p.ANm.resize(p.n);
  p.ALm.resize(p.n);
  p.AFm.resize(p.n);
  p.Jmax.resize(p.n);
  p.Jresp.resize(p.n);
  p.mort.resize(p.n);
  p.mHTL.resize(p.n);

  p.theta.resize(p.n);
  for (int i=0; i<p.n; i++) {
    p.m[i] = _m[i];
    p.ANm[i] = _ANm[i];
    p.ALm[i] = _ALm[i];
    p.AFm[i] = _AFm[i];
    p.Jmax[i] = _Jmax[i];
    p.Jresp[i] = _Jresp[i];
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
  JrespT.resize(p.n);
  //
  //  Init rates:
  //
  rates.JN.resize(p.n);
  rates.JDOC.resize(p.n);
  rates.JL.resize(p.n);
  rates.JF.resize(p.n);
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

extern "C" void printParameters() {
  std::cout << p.n << "\n";
  //std::cout << p.m[0] << "," << p.m[1] << "\n";
  for (int i=0; i<p.n; i++) {
    for (int j=0; j<p.n; j++) {
      std::cout << p.theta[j][i] << ",";
    };
    std::cout << "\n";
  };
};

inline double fTemp(double Q10, double T) {
  return pow(Q10, T/10-1);
};

extern "C" void calcRates(const double& T, const double& N,
			  const double& DOC, const double& L, const double* _B, double* dudt) {
  int i;

  for (i=0; i<p.n; i++)
    B[i] = (0 < _B[i]) ? _B[i] : 0;

  rates.JN[0] = 10;
    //
    // Temperature corrections:
    //
    double f15 = fTemp(1.5, T);
    double f20 = fTemp(2.0, T);
    for (i=0; i<p.n; i++) {
      ANmT[i] = p.ANm[i]*f15;
      JmaxT[i] = p.Jmax[i]*f20;
      JrespT[i] = p.Jresp[i]*f20;
    }
    //
    // Uptakes
    //
    for (i=0; i<p.n; i++) {
      // Uptakes:
      rates.JN[i] = JmaxT[i] * ANmT[i]*N / (JmaxT[i] + ANmT[i]*N*p.rhoCN);
      rates.JDOC[i] = JmaxT[i] * ANmT[i]*DOC / (JmaxT[i] + ANmT[i]*DOC);
      rates.JL[i] = p.epsilonL * JmaxT[i] * p.ALm[i]*L / (JmaxT[i] + p.ALm[i]*L);

      rates.F[i] = 0;
      for (int j = 0; j<p.n; j++) {
	rates.F[i] += p.theta[j][i]*B[j];
      }
      rates.JF[i] = p.epsilonF * JmaxT[i] * p.AFm[i]*rates.F[i] / (JmaxT[i] +p.AFm[i]*rates.F[i]);  
      // Downregulation:
      rates.JNtot[i] = rates.JN[i] + rates.JF[i]/p.rhoCN;
      rates.JLreal[i] = min(rates.JL[i], max(0, rates.JNtot[i]*p.rhoCN - (rates.JF[i]+rates.JDOC[i]-JrespT[i])));
      rates.JCtot[i] = rates.JLreal[i] + rates.JF[i] + rates.JDOC[i] - JrespT[i];
      // Synthesis:
      rates.Jtot[i] = min( rates.JCtot[i], rates.JNtot[i]*p.rhoCN);
      // Losses:
      rates.JCloss_feeding[i] = (1-p.epsilonF)/p.epsilonF*rates.JF[i];
      rates.JCloss_photouptake[i] = (1-p.epsilonL)/p.epsilonL*rates.JLreal[i];
      rates.JNlossLiebig[i] = max(0, rates.JNtot[i]*p.rhoCN-rates.JCtot[i])/p.rhoCN;
      rates.JClossLiebig[i] = max(0, rates.JCtot[i]-rates.JNtot[i]*p.rhoCN); 
      
      rates.JNloss[i] = rates.JCloss_feeding[i]/p.rhoCN +rates. JNlossLiebig[i];
      rates.JCloss[i] = rates.JCloss_feeding[i] +rates. JCloss_photouptake[i] + rates.JClossLiebig[i];
    }
    //
    // Mortality:
    //
    for (i=0; i<p.n; i++) {
      rates.mortpred[i] = 0;
      for (int j = 0; j<p.n; j++) {
	rates.mortpred[i] += p.theta[i][j] * rates.JF[j]*B[j]/(p.epsilonF*p.m[j]*rates.F[j]);
      }
    }
    //
    // Reaction terms:
    //
    dudt[idxN] = 0;
    dudt[idxDOC] = 0;
    for (i=0; i<p.n; i++) {
      double mortloss;
      mortloss = B[i]*(p.mort2*B[i] + p.mHTL[i]);
      dudt[idxN] += (-rates.JN[i]+rates.JNloss[i])*B[i]/p.m[i] + p.remin*mortloss/p.rhoCN;
      dudt[idxDOC] += (-rates.JDOC[i] + rates.JCloss[i])*B[i]/p.m[i] + p.remin*mortloss;
      
      dudt[idxB+i] = (rates.Jtot[i]/p.m[i]  - (p.mort[i] + rates.mortpred[i] + p.mort2*B[i] + p.mHTL[i]))*B[i];
      //Bdt[(B<1e-3) & (dBdt<0)] = 0 # Impose a minimum concentration even if it means loss of mass balance

      //std::cout << rates.JNloss[i] << "," << mortloss << "\n";
    }
};

extern "C" void derivativeChemostat(const double& L, const double& T, const double& d, const double& N0,
				    const double* u, double* dudt) {
  calcRates(T, u[idxN], u[idxDOC ], L, &u[idxB], dudt);
  dudt[idxN] += d*(N0-u[idxN]);
  dudt[idxDOC] += d*(0-u[idxDOC]);
  for (int i=0; i<p.n; i++)
    dudt[idxB+i] += d*(0-u[idxB+i]);
};

int main()
{
  return 0;
};
