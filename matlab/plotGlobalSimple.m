function plotGlobalSimple(sim, iTime)
% Choose the last timestep if none is given:
if (nargin()==1)
    iTime = length(sim.t);
end

% Load grid:
load('../TMs/MITgcm/grid.mat');
x = [x-x(1) ;360]; % adjust x coordinates to map plot
%%
% Prepare timestep:
%
N(:,:,:) = matrixToGrid(sim.N(:,iTime), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
DOC(:,:,:) = double(matrixToGrid(sim.DOC(:,iTime), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
B(:,:,:) = double(matrixToGrid(sum(sim.B(:,:,iTime),2), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));

%%
% Do the plots:

clf
set(gcf,'color','w');


% DOC
%text(0, 1, labels(i),'Units','normalized')
subplot(3,1,1)
panelGlobal(x,y,DOC(:,:,1),'DOC');

% Nitrogen
subplot(3,1,2)
c = panelGlobal(x,y,N(:,:,1),'N');
c.Label.String  = 'Concentration [\mug N l^{-1}]';
box off
%text(0, 1, labels(i),'Units','normalized')

% Plankton
subplot(3,1,3)
panelGlobal(x,y,B(:,:,1),'Plankton');
