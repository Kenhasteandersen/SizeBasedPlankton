p = parametersGlobal;

p.dt = 0.01;
sim = simulateGlobal(p);

p.dt = 0.1;
sim2 = simulateGlobal(p);

p.dt = 0.25;
sim3 = simulateGlobal(p);

%%
load(p.pathGrid);
x = [x-x(1) ;360]; % adjust x coordinates to map plot
B(:,:,:)  = double(matrixToGrid(sum(sim.B(:,:,end),2), [], p.pathBoxes, p.pathGrid));
B2(:,:,:) = double(matrixToGrid(sum(sim2.B(:,:,end),2), [], p.pathBoxes, p.pathGrid));
B3(:,:,:) = double(matrixToGrid(sum(sim3.B(:,:,end),2), [], p.pathBoxes, p.pathGrid));


clf
subplot(3,1,1)
panelGlobal(x,y,log10(B(:,:,1)),'dt=0.01','fast');
caxis([1 4])

subplot(3,1,2)
panelGlobal(x,y,log10(B2(:,:,1)),'dt=0.1','fast');
caxis([1 4])

subplot(3,1,3)
panelGlobal(x,y,log10(B3(:,:,1)),'dt=0.25','fast');
caxis([1 4])
