function plotGlobalWatercolumn(lat,lon,sim,iTime)
if (nargin()==3)
    iTime = length(sim.t);
end

[idx, z] = calcGlobalWatercolumn(lat,lon);


clf
subplot(1,2,1)
loglog(sim.N(idx,iTime)/10,-z, 'b-','linewidth',2)
hold on
loglog(sim.DOC(idx,iTime),-z, 'm-','linewidth',2)
loglog(sum(sim.B(idx,:,iTime),2),-z,'r-','linewidth',2)
loglog(sim.L(idx,iTime)+eps,-z,'y-','linewidth',2)
ylim(-[max(z) min(z)])
xlim([0.01,2000])
legend({'N/10','DOC','B','Light'},'location','northoutside')
ylabel('Depth (m)')
xlabel('Concentration (\mu g)')

subplot(1,2,2)
surface(sim.p.m, -log10(z), squeeze(sim.B(idx,:,iTime)));
set(gca,'xscale','log')
axis tight
shading interp
xlabel('mass (\mu gC)')
colorbar('northoutside')


end