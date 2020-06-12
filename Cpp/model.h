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
/*
 * Globals:
 */
plankton p;
typeRates rates;
type_mass B, ANmT, JmaxT, JFmaxT, JrespT;
const int idxN = 0;
const int idxDOC = 1;
const int idxB = 2;
