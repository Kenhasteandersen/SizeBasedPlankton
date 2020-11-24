function sim = simulateChemostat(p, tEnd)
%
% Load library:
%
loadModel;
%
% Set parameters:
%
setParameters(p)
%
% Simulate:
%
[t, y] = ode23(@derivChemostat, [0 tEnd], [0.1*p.N0, p.DOC0, p.B0], [], p);
%
% Extract result:
%
ixB = 3:p.n+2;
ixTime = find(t>=(max(t)/2),1):length(t);

sim.t = t;
sim.p = p;

sim.Btime = y(:,ixB);
sim.Ntime = y(:,1);
sim.DOCtime = y(:,2);

sim.Bmean = mean(sim.Btime(ixTime,:));
sim.Nmean = mean(sim.Ntime(ixTime,:));
sim.DOCmean = mean(sim.DOCtime(ixTime,:));

end
