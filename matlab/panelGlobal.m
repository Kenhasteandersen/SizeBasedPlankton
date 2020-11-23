%
% Make a panel of a global field given in gridded coordinates:
%
function cbar = panelGlobal(x,y,z, sTitle)

ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','eckert4', ...
    'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
    'off', 'FLineWidth', 2);
ax.XColor = 'white';
ax.YColor = 'white';
axis tight manual
%plabel('PlabelLocation',20, 'PLabelMeridian', 91)
surfacem(y,x ,squeeze(z)');
%shading interp
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
cbar = colorbar('southoutside', 'FontSize',14);
cbar.Label.String  = 'Concentration [\mug C l^{-1}]';
box off
title(sTitle)


