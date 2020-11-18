function gifGlobal(sim)
% Load grid:
load('../TMs/MITgcm/grid.mat');
xp = [x-x(1) ;360]; % adjust x coordinates to map plot

labels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
%% 
% Divide into pico, nano and micro (to be done more rigorously)

pB = squeeze(sum(sim.B(:,1:3,:),2)); %pico
nB = squeeze(sum(sim.B(:,4:7,:),2)); %nano
mB = squeeze(sum(sim.B(:,8:10,:),2)); %micro

% taking last 12 months if sim.t >12
pB = pB(:,end-11:end);
nB = nB(:,end-11:end);
mB = mB(:,end-11:end);

% Prepare all timesteps:
%
for i = 1:length(sim.t)
    Nplot(:,:,:,i) = matrixToGrid(sim.N(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    DOCplot(:,:,:,i) = double(matrixToGrid(sim.DOC(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    pBplot(:,:,:,i) = double(matrixToGrid(pB(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    nBplot(:,:,:,i) = double(matrixToGrid(nB(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    mBplot(:,:,:,i) = double(matrixToGrid(mB(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
end

%%
% Do the plots:

% DOC
DOClim = [0 max(DOCplot(:,:,1,:),[],'all')];

 h = figure('Position', [300, 300, 700, 400]);
    set(gcf,'color','w');    
for i = 1:12
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
         'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
         'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(DOCplot(:,:,1,i))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = '[\mug C l^{-1}]';
    caxis(DOClim);
    box off
    title('DOC')
    text(0, 1, labels(i),'Units','normalized')
    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'DOC.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'DOC.gif','gif','WriteMode','append');
    end
    clf
end

close(h)
% Nitrogen
Nlim = [0 max(Nplot(:,:,1,:),[],'all')];

 h = figure('Position', [400, 400, 700, 400]);
    set(gcf,'color','w');  
for i = 1:12
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
         'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
         'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(Nplot(:,:,1,i))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    c = colorbar('southoutside', 'FontSize',14);
    c.Label.String  = '[\mug N l^{-1}]';
    caxis(Nlim);
    box off
    title('N')
    text(0, 1, labels(i),'Units','normalized')
    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'N.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'N.gif','gif','WriteMode','append');
    end
    clf
end

close(h)
% Plankton
pBlim = [0 max(pBplot(:,:,1,:),[],'all')];
nBlim = [0 max(nBplot(:,:,1,:),[],'all')];
mBlim = [0 max(mBplot(:,:,1,:),[],'all')];

h = figure('Position', [0, 0, 1400, 400]); 
set(gcf,'color','w');    

for i=1:12
   subplot(1,3,1)
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual

    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(pBplot(:,:,1,i))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    %if i==10
        %caxis([0 16])
        c = colorbar('southoutside', 'FontSize',14);
        c.Label.String  = '[\mu g C l^{-1}]';
        caxis(pBlim);
    %end
    box off
    title('Pico')
    text(0, 1, labels(i),'Units','normalized')
    
    
    subplot(1,3,2)
    set(gcf,'color','w');
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(nBplot(:,:,1,i))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
    %if i==10
        c = colorbar('southoutside', 'FontSize',14);
        c.Label.String  = '[\mu g C l^{-1}]';
        caxis(nBlim);
    %end
    box off
    title('Nano')
    %text(0, 1, labels(i),'Units','normalized')
    
    subplot(1,3,3)
    set(gcf,'color','w');
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', 'Grid', 'on', 'Frame', 'on',...
        'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,xp ,squeeze(mBplot(:,:,1,i))');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
        c = colorbar('southoutside', 'FontSize',14);
        c.Label.String  = '[\mu g C l^{-1}]';
        caxis(mBlim);
    box off
    title('Micro')
    %text(0, 1, labels(i),'Units','normalized')
    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,'Biomass_PNM.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'Biomass_PNM.gif','gif','WriteMode','append');
    end
    clf
end
end
