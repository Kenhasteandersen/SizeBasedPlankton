function plotGlobal(sim, iTime)
% Choose the last timestep if none is given:
if (nargin()==1)
    iTime = length(sim.t);
end

n = double(sim.p.n);
% Load grid:
load('../TMs/MITgcm/grid.mat');

%%
% Prepare all timesteps:
%
for i = 1:length(sim.t)
    Nplot(:,:,:,i) = matrixToGrid(sim.N(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    DOCplot(:,:,:,i) = double(matrixToGrid(sim.DOC(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    for j = 1:n
        Bplot(:,:,:,j,i) = double(matrixToGrid(sim.B(:,j,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    end
end

%%
% Do the plots:
%
figure
for i=1:n
    subplot(2,floor(n/2),i)
    surface(x,y,Bplot(:,:,1,i,iTime)')
    title(['B',num2str(i)])
    c = colorbar;
    %caxis([0 1.5]);
    shading flat
end

figure
surface(x,y,DOCplot(:,:,1,iTime)')
title('DOC')
c=colorbar;
shading flat

figure
surface(x,y,Nplot(:,:,1,iTime)')
title('N')
c=colorbar;
shading flat
