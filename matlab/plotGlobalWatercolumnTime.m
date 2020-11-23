function plotGlobalWatercolumnTime(lat,lon,sim)

[idx, z] = calcGlobalWatercolumn(lat,lon);
t = sim.t;


clf
subplot(3,1,1)
surface(t,-log10(z), log10(double(sim.N(idx,:))))
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('log10 Depth (m)')
set(gca,'yscale','log')
shading interp
axis tight
colorbar

subplot(3,1,2)
surface(t,-log10(z), log10(double(sim.DOC(idx,:))))
title(['DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('log10 Depth (m)')
%xlabel('Concentration (\mugC l^{-1})')
set(gca,'yscale','log')
shading interp
axis tight
colorbar

subplot(3,1,3) 
surface(t,-log10(z), squeeze(log10(sum(double(sim.B(idx,:,:)),2))))
title(['Plankton, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('log10 Depth (m)')
xlabel('Time (days)')
set(gca,'yscale','log')
shading interp
axis tight
colorbar
caxis([-1 2])

end