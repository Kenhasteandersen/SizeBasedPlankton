function plotSpectrum(sim)

m = sim.p.m;



loglog(m, sim.Bmean, 'k', 'linewidth',3)
xlim([min(m), max(m)])
ylim([1, 500])

xlabel('Carbon mass (\mu gC)')
ylabel('Biomass (\mu gC/l)')