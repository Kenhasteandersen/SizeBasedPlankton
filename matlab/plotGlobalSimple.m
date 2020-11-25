function plotGlobalSimple(sim, iTime)
% Choose the last timestep if none is given:
if (nargin()==1)
    iTime = length(sim.t);
end

sType = 'fast';
%
% Load grid:
%
load(sim.p.pathGrid);
%
% Prepare timestep:
%
N(:,:,:) = matrixToGrid(sim.N(:,iTime), [], sim.p.pathBoxes, sim.p.pathGrid);
DOC(:,:,:) = double(matrixToGrid(sim.DOC(:,iTime), [], sim.p.pathBoxes, sim.p.pathGrid));
B(:,:,:) = double(matrixToGrid(sum(sim.B(:,:,iTime),2), [], sim.p.pathBoxes, sim.p.pathGrid));
C=calcGlobalCnet(sim,iTime);
Cnet = double(matrixToGrid(C, [], sim.p.pathBoxes, sim.p.pathGrid));


%
% Do the plots:
%
clf
set(gcf,'color','w');

% DOC
%text(0, 1, labels(i),'Units','normalized')
subplot(4,1,1)
panelGlobal(x,y,DOC(:,:,1),'DOC',sType);

% Nitrogen
subplot(4,1,2)
c = panelGlobal(x,y,N(:,:,1),'N',sType);
c.Label.String  = 'Concentration [\mug N l^{-1}]';

% Plankton
subplot(4,1,3)
panelGlobal(x,y,log10(B(:,:,1)),'Plankton (log10)',sType);
caxis([1 3])

subplot(4,1,4)
panelGlobal(x,y,log10(Cnet(:,:,1)),'Net primary production TOP LAYER ONLY (log10)', sType);
caxis([0 2])