%
% Calculated primary production
%
function Cnet = calcGlobalCnet(sim, iTime)

Cnet = zeros(size(sim.N,1),1);
for k = 1:size(sim.N,1)
    Cnet(k) = calcCnet(sim.p, sim.N(k,iTime), sim.DOC(k,iTime), sim.B(k,:,iTime), sim.L(k,iTime), sim.T(k,iTime));
end

