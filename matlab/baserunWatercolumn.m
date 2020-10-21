function sim = baserunWatercolumn
%%
p = parametersWatercolumn;
sim = simulateWatercolumn(p, true, p.diff*ones(1,p.nGrid), 110);
plotWatercolumn(sim)


