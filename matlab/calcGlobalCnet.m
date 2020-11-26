%
% Calculated primary production
%
function Cnet = calcGlobalCnet(sim, iTime)

Cnet = 0*squeeze(sim.N(:,:,:,1));
for i = 1:size(sim.x)
    for j = 1:size(sim.y)
        for k = 1:size(sim.z)
            Cnet(i,j,k) = calcCnet(sim.p, sim.N(i,j,k,iTime), ...
                sim.DOC(i,j,k,iTime), squeeze(sim.B(i,j,k,:,iTime))', sim.L(i,j,k,iTime), sim.T(i,j,k,iTime));
        end
    end
end

