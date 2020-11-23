function plotGlobal(sim, iTime)
% Choose the last timestep if none is given:
if (nargin()==1)
    iTime = length(sim.t);
end

n = double(sim.p.n);
% Load grid:
load('../TMs/MITgcm/grid.mat');
xp = [x-x(1) ;360]; % adjust x coordinates to map plot
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

% DOC

 h = figure('Position', [300, 300, 700, 400]);
    set(gcf,'color','w');    
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
         'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
         'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(DOCplot(:,:,1,iTime))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = 'Concentration [\mug C l^{-1}]';
    box off
    title('DOC')
    %text(0, 1, labels(i),'Units','normalized')


% Nitrogen
 h = figure('Position', [400, 400, 700, 400]);
    set(gcf,'color','w');    
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
         'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
         'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(Nplot(:,:,1,iTime))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = 'Concentration [\mug N l^{-1}]';
    box off
    title('N')
    %text(0, 1, labels(i),'Units','normalized')

% Plankton
for ii = 1:n
    figure('Position', [500, 500, 700, 400]);
    set(gcf,'color','w');    
     ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
         'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
         'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(Bplot(:,:,1,ii,iTime))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    %if i==10
        %caxis([0 12])
        c = colorbar('southoutside', 'FontSize',14);
        c.Label.String  = 'Concentration [\mug C l^{-1}]';
    %end
    box off
    title(['B ',num2str(ii)])
    %text(0, 1, labels(i),'Units','normalized')
end
