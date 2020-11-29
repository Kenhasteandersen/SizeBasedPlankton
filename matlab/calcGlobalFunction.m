function sim = calcGlobalFunction(sim)
%
% Primary production:
%
for i = 1:length(sim.t)
    sim.Cnet(:,:,:,i) = calcGlobalCnet(sim,i);
end
%
% Global totals
%
    function tot = calcTotal(u)
        tmp = squeeze(u(:,:,:,1));
        tot = sum(u(ix).*dv(ix));
    end

% Get grid volumes:
load(sim.p.pathGrid,'dv');

ix = ~isnan(sim.N(:,:,:,1)); % Find all relevant grid cells
ix = ix(:);

for i = 1:length(sim.t)
    sim.Ntotal(i) = calcTotal(sim.N(:,:,:,i));
    sim.DOCtotal(i) = calcTotal(sim.DOC(:,:,:,i));
    sim.Btotal(i) = 0;
    for j = 1:sim.p.n
        sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(:,:,:,j,i));
    end
    sim.CnetTotal(i) = calcTotal(sim.Cnet(:,:,:,i));
end

end