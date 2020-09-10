function [t, y] = simulateChemostat(p)
%
% Load library:
%
if libisloaded("model")
    unloadlibrary("model")
end
loadlibrary('../Cpp/model.so','../Cpp/model.h')
%
% Set parameters:
%
calllib('model','setParameters', ...
    p.n, p.m, p.rhoCN, p.epsilonL, p.epsilonF, ...
    p.ANm, p.ALm, p.AFm, ...
    p.Jmax, p.JFmaxm, p.Jresp, p.Jloss_passive_m, ...
    p.theta, p.mort, p.mort2, p.mortHTL*p.mortHTLm, ...
    p.remin, p.remin2, p.cLeakage)
%
% Simulate:
%
[t, y] = ode45(@derivChemostat, [0 1000], [p.N0, p.DOC0, p.B0], [], p);

end
