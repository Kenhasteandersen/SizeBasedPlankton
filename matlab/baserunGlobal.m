p = parametersGlobal();

load('../TMs/globalInit.mat'); % Load decent initial conditions
sim = simulateGlobal(p,sim);

close all
figure
plotGlobalSimple(sim);

figure
plotGlobalWatercolumnTime(60,0,sim);

%plotGlobal(sim);
%gifGlobal(sim);