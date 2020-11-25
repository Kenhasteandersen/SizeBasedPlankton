function plotGlobalWatercolumn(lat,lon,sim,iTime)
if (nargin()==3)
    iTime = length(sim.t);
end

[idx, z] = calcGlobalWatercolumn(lat,lon, sim);


clf
subplot(1,2,1)
plot(sim.N(idx,iTime)/10,-z, 'b-','linewidth',2)
hold on
plot(sim.DOC(idx,iTime),-z, 'm-','linewidth',2)

for i=1:length(idx)
    Bpnm(i,:) = calcPicoNanoMicro(sim.B(idx(i),:,iTime), sim.p);
end
for i = 1:3
    loglog(Bpnm(:,i),-z,'r-','linewidth',i)
end
plot(sim.L(idx,iTime)+eps,-z,'y-','linewidth',2)
ylim(-[max(z) 0])
xlim([0.01,2000])
legend({'N/10','DOC','Bpico','Bnano','Bmicro','Light'},'location','northoutside')
ylabel('Depth (m)')
xlabel('Concentration ({\mu}g/L)')
set(gca,'xscale','log')


subplot(1,2,2)
surface(sim.p.m, -z, squeeze(sim.B(idx,:,iTime)));
set(gca,'xscale','log')
axis tight
shading interp
xlabel('mass (\mu gC)')
colorbar('northoutside')


end