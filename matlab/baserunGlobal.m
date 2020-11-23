%
% Make a basic run of the global model
%
p = parametersGlobal();

load('../TMs/globalInit.mat'); % Load decent initial conditions
p.tSave = 10; % Save every 10th day
sim = simulateGlobal(p,sim); % Simulate
%
% Plots:
%
close all
figure
plotGlobalSimple(sim);

figure
plotGlobalWatercolumnTime(60,10,sim);
%
% CPU-heavy plots:
%
%plotGlobal(sim);
%gifGlobal(sim);