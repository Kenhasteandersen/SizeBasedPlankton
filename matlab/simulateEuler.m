%
% Simulate one system with Euler-integration:
%
function sim = simulateEuler(p, y, dt, tEnd)
%%
% Load library:
%
loadModel;
%
% Set parameters:
%
setParameters(p);
%%
% Simulate:
%

dydt = 0*y;
y = calllib('model','simulateEuler', y, dydt, ...
    double(p.L), double(p.T), double(dt), double(tEnd));

%%
% Extract result:
%
sim.t = tEnd;
sim.p = p;

sim.Ntime = y(1);
sim.DOCtime = y(2);
sim.Btime = y(3:(2+p.n));

%sim.Bmean = mean(sim.Btime(ixTime,:));
%sim.Nmean = mean(sim.Ntime(ixTime,:));
%sim.DOCmean = mean(sim.DOCtime(ixTime,:));

end
