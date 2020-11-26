function plotGlobalWatercolumn(lat,lon,sim,iTime)
if (nargin()==3)
    iTime = length(sim.t);
end

idx = calcGlobalWatercolumn(lat,lon, sim);
z = -sim.z(idx.z);

if length(z)>0
    clf
    subplot(1,2,1)
    plot(squeeze(sim.N(idx.x, idx.y, idx.z, iTime))/10, z, 'b-','linewidth',2)
    hold on
    plot(squeeze(sim.DOC(idx.x, idx.y, idx.z,iTime)), z, 'm-','linewidth',2)
    
    for i=1:length(idx.z)
        Bpnm(i,:) = calcPicoNanoMicro(squeeze(sim.B(idx.x, idx.y, idx.z(i),:,iTime)), sim.p);
    end
    for i = 1:3
        loglog(Bpnm(:,i), z,'r-','linewidth',i)
    end
    plot(squeeze(sim.L(idx.x, idx.y, idx.z,iTime))+eps, z,'y-','linewidth',2)
    ylim([min(z) 0])
    xlim([0.01,2000])
    legend({'N/10','DOC','Bpico','Bnano','Bmicro','Light'},'location','northoutside')
    ylabel('Depth (m)')
    xlabel('Concentration ({\mu}g/L)')
    set(gca,'xscale','log')
    
    
    subplot(1,2,2)
    surface(sim.p.m,  z, squeeze(sim.B(idx.x, idx.y, idx.z,:,iTime)));
    set(gca,'xscale','log')
    axis tight
    shading interp
    xlabel('mass (\mu gC)')
    colorbar('northoutside')
end

end