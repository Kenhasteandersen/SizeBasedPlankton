p = parametersGlobal();

load('../TMs/globalInit.mat'); % Load decent initial conditions
sim = simulateGlobal(p, sim);

close all
plotGlobal(sim);
