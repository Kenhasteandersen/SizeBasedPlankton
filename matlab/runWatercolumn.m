%
% Run water column with gotm annual simulation
%

p = parameters;

p.nGrid = 30;
p.diff = 1;
%%

p.tEnd = 200;
p.dt = 0.05;
p.Lamplitude = 0;
%[p,sol] = modelWatercolumn(p,0.01);

%p.epsilonL = 1;
[p,sol] = modelWatercolumnImplicit(p);
%[p,sol] = modelWatercolumnImplicitPC(p);
