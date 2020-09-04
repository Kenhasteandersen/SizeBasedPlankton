void setParameters(
    const int* _n, 
    const double* _m,
    const double* _rhoCN, 
    const double* _epsilonL, 
    const double* _epsilonF,
    const double* _ANm, 
    const double* _ALm, 
    const double* _AFm,
    const double* _Jmax, 
    const double* _JFmax,
    const double* _Jresp,
    const double* _Jloss_passive,
    const double* _theta,
    const double* _mort,
    const double* _mort2,
    const double* _mortHTL,
    const double* _remin,
    const double* _remin2,
    const double* _cLeakage
);

//void calcRates(const double* T, const double* L, const double* u, double* dudt);

void derivativeChemostat(const double* L, const double* T, const double* d, const double* N0,
                                   const double* u, double* dudt);
/*
void simulateWaterColumnFixed(const double& L0, const double& T, 
                                        const double* Diff, const double& N0,
                                        const double& tEnd, const double& dt,
                                        const int& nGrid, const double* x, double* u);
 
 */