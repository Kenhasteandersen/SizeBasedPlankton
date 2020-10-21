function plotWatercolumn(sim)

%%
figure(1)
clf
surface(sim.p.m, sim.x, squeeze(sim.B(:,:,end))')
set(gca, 'xscale','log')
xlabel('Mass (\mu gC)')
ylabel('Depth (m)')
axis tight
shading interp

figure(2)
clf
subplot(5,1,1)
surface(sim.t, sim.x, sim.N)
shading interp
axis tight
colorbar

subplot(5,1,2)
surface(sim.t, sim.x, sim.DOC)
shading interp
axis tight
colorbar

subplot(5,1,3)
surface(sim.t, sim.x, sim.Bpico)
shading interp
axis tight
colorbar

subplot(5,1,4)
surface(sim.t, sim.x, sim.Bnano)
shading interp
axis tight
colorbar

subplot(5,1,5)
surface(sim.t, sim.x, sim.Bmicro)
shading interp
axis tight
colorbar

figure(3)
clf
plot(sim.Bpico(:,end), sim.x,'k-')
hold on
plot(sim.Bnano(:,end), sim.x, 'k-', 'linewidth',2)
plot(sim.Bmicro(:,end), sim.x, 'k-', 'linewidth',3)

plot(sim.N(:,end)/10, sim.x, 'b-','linewidth',2)
plot(sim.DOC(:,end)*10, sim.x, 'm-', 'linewidth',2)

ylim([-sim.p.depth,0])
ylabel('depth (m)')
xlabel('(\mu gC/l)')

legend('Pico','Nano','Micro','N/10','DOC*10', 'location','southeast')

