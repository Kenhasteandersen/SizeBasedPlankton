function plotGlobalSimple(sim, iTime)
% Choose the last timestep if none is given:
if (nargin()==1)
    iTime = length(sim.t);
end

sType = 'fast';
%
% Calc primary production:
%
Cnet = calcGlobalCnet(sim,iTime);


%
% Do the plots:
%
clf
set(gcf,'color','w');

% DOC
%text(0, 1, labels(i),'Units','normalized')
subplot(4,1,1)
panelGlobal(sim.x,sim.y,sim.DOC(:,:,1,iTime),'DOC',sType);

% Nitrogen
subplot(4,1,2)
c = panelGlobal(sim.x,sim.y,sim.N(:,:,1,iTime),'N',sType);
c.Label.String  = 'Concentration [\mug N l^{-1}]';

% Plankton
subplot(4,1,3)
panelGlobal(sim.x,sim.y,log10(sim.B(:,:,1,iTime)),'Plankton (log10)',sType);
caxis([1 3])

subplot(4,1,4)
panelGlobal(sim.x, sim.y, log10(Cnet(:,:,1)),'Net primary production TOP LAYER ONLY (log10)', sType);
caxis([0 2])