%
% Compare Euler and ode45 integation:
%
function compareODEvsEuler(p, dt)
%p = parametersChemostat();
p.d = 0;

tEnd = 365;

tic
simODE45 = simulateChemostat(p, tEnd);
toc

tic
simEuler = simulateEuler(p, [p.N0, p.DOC0, p.B0], dt, tEnd);
toc

clf
loglog(p.m, simODE45.Btime(end,:))
hold on
loglog(p.m, simEuler.Btime,'r--')

