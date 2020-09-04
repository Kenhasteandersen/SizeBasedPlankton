function p = parameters()
%
% Numerical parameters:
%
p.nGrid = 30; % No. of grid points
p.dt = 0.1; % dt for fixed time-step algorithms
p.tEnd = 365;
%
% Define cell groups:
%
p.n = 10; % No of groups
p.m = logspace(-8,1,p.n)';  % Mass bins in mugC

p.rhoCN = 5.68; % C:N mass ratio
p.epsilonL = 0.8; % Light uptake efficiency
p.epsilonF = 0.8; % Assimilation efficiency
%
% Cell wall fraction of mass:
%
p.c = 0.0015; % the constant is increased a bit to limit the lower cell size
nu = p.c * p.m.^(-1/3);
%
% Clearance rates:
%
p.AN = 0.00004; % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g).^(-1/3) / 1.5 (g/cm); Andersen et al 2015
p.AL = 0.000914; % if using Andys shading formula for non-diatoms
p.cL = 21; % if using Andys shading formula for non-diatoms
p.AF = 0.018;  %  Fits to TK data for protists

p.ANm = p.AN*p.m.^(1/3);
p.ALm = p.AL*p.m.^(2/3) .* (1-exp(- p.cL*p.m.^(1/3) ));  % shading formula
p.AFm = p.AF*p.m;
%
% Prey encounter
%
p.beta = 500;
p.sigma = 1.3;
p.theta = zeros(p.n, p.n);
for i = 1:p.n
    p.theta(i,:) = exp( -(log(p.m(i)./p.m/p.beta)).^2/(2*p.sigma^2) );
end
%
% Metabolism:
%
p.alphaJ = 1.5;
p.Jmax = p.alphaJ * p.m .* (1-nu); % mugC/day
p.cR = 0.1;
p.Jresp = p.cR*p.alphaJ*p.m;
%
% Losses:
%
p.mort = 0*0.005*(p.Jmax./p.m) .* p.m.^(-1/4);
p.mort2 = 0.0002*p.n;
p.mortHTL = 0.2;
p.mHTL = max(p.m)/p.beta; % Bins affected by HTL mortality

p.remin = 0.0; % fraction of mortality losses reminerilized to N and DOC
%
% Biogeochemical model:
%
p.T = 10;   % Temperature
p.N0 = 150; % Deep nutrient levels
p.Depth = 100;
%
% Initial conditions:
%
p.DOC0 = 0;
p.B0 = 10*ones(1, p.n);
%
% Light:
%
p.L = 300;  % PAR, mu E/m2/s
p.Lamplitude = 0.8;

p.diff = 50; % Diffusion for chemostat

p.tEnd = 365; % Simulation length (days)
