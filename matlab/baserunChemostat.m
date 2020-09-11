function baserunChemostat

p = parametersChemostat;
sim = simulateChemostat(p);
plotSpectrum(sim)

end