%
% Compare Euler and ode45 integation:
%
function [simODE45 simEuler] = compareODEvsEuler(p, dt)
%p = parametersChemostat();
%unloadlibrary("model")
p.d = 0.0;

tEnd = 365;

tic
simODE45 = simulateChemostat(p, tEnd);
toc

tic
simEuler = simulateEuler(p, [0.1*p.N0, p.DOC0, p.B0], dt, tEnd);
toc

%u = [p.N0, p.DOC0, p.B0];
%for i = 1:365
%    simEuler = simulateEuler(p, u, dt, 1);
%    u = [simEuler.Ntime simEuler.DOCtime simEuler.Btime];
%end

clf
loglog(p.m, simODE45.Btime(end,:))
hold on
loglog(p.m, simEuler.Btime,'r--')

