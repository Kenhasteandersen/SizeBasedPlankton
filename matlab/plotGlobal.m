load('../Data/global_results/euler_global_0yr_n4.mat') 
load('../TMs/MITgcm/grid.mat');

%%

for i = 1:12;
    Nplot(:,:,:,i) = matrixToGrid(Nm(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    DOCplot(:,:,:,i) = double(matrixToGrid(DOCm(:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    for j = 1:p.n
        Bplot(:,:,:,j,i) = double(matrixToGrid(Bmatm(:,j,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'));
    end

end


%%
m = 12; % month
figure
for i=1:4
    subplot(2,2,i)
    surface(x,y,Bplot(:,:,1,i,m)')
    title(['B',num2str(i)])
    c = colorbar;
    %caxis([0 1.5]);
    shading flat
end

figure
surface(x,y,DOCplot(:,:,1,m)')
title('DOC')
c=colorbar;
shading flat

figure
surface(x,y,Nplot(:,:,1,m)')
title('N')
c=colorbar;
shading flat
