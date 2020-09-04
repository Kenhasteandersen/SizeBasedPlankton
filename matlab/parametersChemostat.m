
function p = parametersChemostat()
%
% Define parameters:
%
p.n = int32(25);
p.m = logspace(-8.5,1,p.n);

p.rhoCN = 5.68; % C:N mass ratio
p.epsilonL = 0.9; % Light uptake efficiency
p.epsilonF = 0.8; % Assimilation efficiency
p.cLeakage = 0.00015; % passive leakage of C and N
%
% Cell wall fraction of mass:
%
p.c = 0.0015; % the constant is increased a bit to limit the lower cell size
nu = p.c * p.m.^(-1/3);
%
% Clearance rates:
%
%factor = (1e-6)^(1/3)/1.5
p.AN = 0.00012; %0.00004 % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
p.cN = 0.1;
p.AL = 0.000914; % if using Andys shading formula for non-diatoms
p.cL = 21; % if using Andys shading formula for non-diatoms
p.AF = 0.018;  %  Fits to TK data for protists
p.cF = 0.6; % Just a guess
%
% Calc rates as a function of m:
%
p.ANm = p.AN*p.m.^(1/3) ./ (1 + p.cN*p.m.^(-1/3));
p.ALm = p.AL*p.m.^(2/3) .* (1-exp(- p.cL*p.m.^(1/3) ));  % shading formula
p.AFm = p.AF*p.m;
p.Jloss_passive_m = p.cLeakage * p.m.^(2/3); % in units of C
p.JFmaxm = p.cF*p.m.^(2/3);
%
% Prey encounter
%
p.beta = 500;
p.sigma = 1.3;
p.theta = zeros(p.n, p.n);
for i = 1:p.n
    %p.theta[i,] = phi(p.m[i]/p.m, p.beta, p.sigma)   % (predator, prey)
    for j = 1:p.n
        p.theta(i,j) = exp( -(log(p.m(i)/p.m(j)/p.beta))^2/(2*p.sigma^2));
    end
end

%
% Metabolism:
%
p.alphaJ = 1.5; % per day
p.Jmax = p.alphaJ * p.m .* (1-nu); % mugC/day
p.cR = 0.1;
p.Jresp = p.cR*p.alphaJ*p.m;
%
% Losses:
%
p.mort = 0*0.005*(p.Jmax./p.m) .* p.m.^(-1/4);
p.mort2 = 0.0002*p.n;
p.mortHTL = 0.1;
p.mHTL = max(p.m)/p.beta^1.5; % Bins affected by HTL mortality
p.mortHTLm = 1./(1+(p.m./p.mHTL).^(-2));

p.remin = 0.0; % fraction of mortality losses reminerilized to N and DOC
p.remin2 = 1; % fraction of virulisus remineralized to N and DOC

p.T = 10;
p.latitude=0;
p.L = 100;
%
% Initial conditions:
%
p.N0 = 150;
p.DOC0 = 0;
p.B0 = ones(1,p.n);

p.d = 0.05;

%
% Prey size function integrated over size groups:
%
    function fTot = Phi(z, p, Delta)
        m = logspace(-12,3,1000);
        dm = diff(m);
        m = m(2:1000);
        
        fPrey = function(m, w0, Delta) {
            integrate( function(logw) phi(m/exp(logw),beta=p$beta,sigma=p$sigma), log(w0/sqrt(Delta)), log(w0*sqrt(Delta)))$value
                }
                
                fTot = function(m0,w0, Delta) {
                    ix = m>m0/sqrt(Delta) & m<m0*sqrt(Delta)
                    sum(
                    as.numeric(
                    lapply( m[ix], function(m) fPrey(m,w0,Delta)/m)
                        )*dm[ix] / log(Delta)^2
                        )
                        }
                        
                        return(fTot(1,1/z,Delta))
                        }
                        
                        
                        
                    end
