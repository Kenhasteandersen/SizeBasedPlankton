%% Global run of size based unicelluar plankton model
% Transport Matrix used: MITgcm 2.8 degrees
% TM must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and be put into the location '../TMs'
%
% ../Cpp/model.cpp must be compiled locally to create model.so


clear all
close all
%% Load Initial January TM and configuations
load('../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_01.mat');
load('../TMs/MITgcm/grid.mat');
load('../TMs/MITgcm/config_data.mat');
load('../TMs/MITgcm/Matrix5/Data/boxes.mat')

% Load paths for switching TM
loadPath = '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
loadPath1 =  '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_';

% Preparing for timestepping. Using 43200 sec timesteps (0.5 days)
Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);

% number of elements in each layer

lyr_n = permute(sum(sum(bathy)),[3,1,2]);

lyr_ind = [1 ; cumsum(lyr_n(1:end-1))+1];

lyr_end = cumsum(lyr_n);

%% Parameter setup and initial conditions
n = 4;  %number of size classes plankton

p = parametersGlobal(n);

% Load library:
%
if libisloaded("model")
    unloadlibrary("model")
end
loadlibrary('../Cpp/model.so','../Cpp/model.h')
%
% Set parameters:
%
calllib('model','setParameters', ...
    p.n, p.m, p.rhoCN, p.epsilonL, p.epsilonF, ...
    p.ANm, p.ALm, p.AFm, ...
    p.Jmax, p.JFmaxm, p.Jresp, p.Jloss_passive_m, ...
    p.theta, p.mort, p.mort2, p.mortHTL*p.mortHTLm, ...
    p.remin, p.remin2, p.cLeakage)

p.nb = nb;

simtime = 730*1; %simulation time in half days

month = 0;
YEAR = 0;

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
mon2 = [31 28 31 30 31 30 31 31 30 31 30 31]*2;
ticks = 1+2*cumsum(mon);

%if no initial conditions available:

%N = zeros(nx,ny,nz) + p.N0;
DOC = zeros(nx,ny,nz) + p.DOC0;
 
% Load initial cond.
load('../Data/N_oceandata.mat')

 
B = zeros(nx,ny,nz,p.n)+0.1; %biomass
 
for i = 1:p.n
    
    B(:,:,1:5,i) = B(:,:,1:5,i)+p.B0(i);
    
end


% Convert from grid form to matrix form
N   = gridToMatrix(N, [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
DOC = gridToMatrix(DOC, [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat'); 

Bmat = zeros(nb,p.n);
for i = 1:p.n
    Bmat(:,i) = gridToMatrix(B(:,:,:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
end

% save matrices half daily 
Nd = single(zeros(nb,simtime));
DOCd = single(zeros(nb,simtime));

Bmatd= single(zeros(nb,p.n,simtime));

% save matrices monthly
Nm = single(zeros(nb,12));
DOCm = single(zeros(nb,12));

Bmatm= single(zeros(nb,p.n,12));







%% Load environmental variables 
%choose which temperature to use.

load('../TMs/MITgcm/BiogeochemData/Theta_bc.mat')
% GCM
% for i = 1:12
% TmatGCM(:,i) = gridToMatrix(Tgcm(:,:,:,i), [], '../bin/MITgcm_ECCO/Matrix1/Data/boxes.mat', '../bin/MITgcm_ECCO/grid.mat');
% end

%Biochemical temperature
Tmat = zeros(nb,12);
for i = 1:12
Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
end


%% Run TM
elapsed_time = zeros(1,simtime);
tic
for i=1:simtime
    telapsed = tic;


    
    % Test for time to change monthly TM
    if ismember(i, ticks+730*YEAR)% 
        month = month + 1;
        TIMESTEP = 0;
        % Reset to Jan
        if month > 12
            month = 1;
        end
        % Load TM
        if month < 10
            load(strcat(loadPath, num2str(month), '.mat'));
            disp(strcat(loadPath, num2str(month), '.mat'))
        else
            load(strcat(loadPath1, num2str(month), '.mat'));
            disp(strcat(loadPath1, num2str(month), '.mat'))

        end
        
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
        
        % Set monthly mean temperature
        p.T = Tmat(:,month);  
        
     
    end
    
    p.i = i/2;                                        %set time of year [days] (for light)
    
    for j=1:length(lyr_n)    
        p.LYR = lyr_ind(j):lyr_end(j);        %current layer indices
        p.phi = Ybox(p.LYR);              %set latitude
        p.I0 = daily_insolation(0,p.phi,p.i,1);  %set light
        p.L(p.LYR) =  p.EinConv*p.PARfrac*p.I0*exp(-p.kw*z(j)); % Light w. latitude
    end
    
    p.L(p.L<1)=0;
    
    
             
       
      
%     tODE = tic;   
    [t, u] = ode23(@derivGlobal, [0:0.25:0.5], [N; DOC; Bmat(:)], [], p);

%disp('ode finished')
%       toc(tODE)
    % only final results for now....
    N      = u(end,1:nb)';
    DOC    = u(end,nb+1:2*nb)';
    Bmat   = reshape(u(end,2*nb+1:end),[],p.n);
    
    % Remove small concentrations
    N(N<1E-6) = 0;
    DOC(DOC<1E-6) = 0;
    Bmat(Bmat<1E-6) = 0;
        
 
    TIMESTEP = TIMESTEP+1;
    % Transport
    N   = Aimp * ( Aexp * N);
    DOC = Aimp * ( Aexp  * DOC);

    for k = 1:size(Bmat,2)   
        Bmat(:,k) =  Aimp * ( Aexp  * Bmat(:,k));
    end
    
    elapsed_time(i) = toc(telapsed);
     %save timeseries
% Daily:     (too large)
%     Nd(:,i) = N;
%     DOCd(:,i) = DOC;
%     Bmatd(:,:,i)= Bmat;
    
% Monthly (mean):
     Nm(:,month) = (Nm(:,month)*(TIMESTEP-1) + single(N))/TIMESTEP;
     DOCm(:,month) = (DOCm(:,month)*(TIMESTEP-1) + single(DOC))/TIMESTEP;
     Bmatm(:,:,month)= (Bmatm(:,:,month)*(TIMESTEP-1) + single(Bmat))/TIMESTEP;
     
     
     
     if ismember(i,cumsum(mon2)+730*YEAR)
         
         save(['../Data/global_results/global_',num2str(YEAR),'yr_n',num2str(n),'.mat'],'Bmatm','Nm','DOCm','elapsed_time','p')
         
     end

    if mod(i,30) == 0
       disp(['t=',num2str(i)]);
       
    end
     if mod(i,730)==0
        
        save(['../Data/global_results/global_',num2str(YEAR),'yr_n',num2str(n),'.mat'],'Bmatm','Nm','DOCm','elapsed_time','p')
        %reset monthly results for new year
        Nm = single(zeros(nb,12));
        DOCm = single(zeros(nb,12));

        Bmatm= single(zeros(nb,p.n,12));
        
        YEAR = YEAR + 1;        
    end
end
solvingtime = toc;
disp(['solvingtime = ',num2str(solvingtime/3600), ' hours']);

%% converting to grid form for plotting (december)
% N = double(matrixToGrid(Nm(:,12), [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat'));
% DOC = double(matrixToGrid(DOCm(:,12), [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat'));
% 
% for i = 1:p.n
%     B(:,:,:,i) = double(matrixToGrid(Bmatm(:,i,12), [], '../../bin/MITgcm/Matrix5/Data/boxes.mat', '../../bin/MITgcm/grid.mat'));
% end
% 
% %%
% figure
% for i = 1:double(p.n)
%     subplot(3,2,i)
%     surface(x,y,B(:,:,1,i)')
%     title(['Biomass ',num2str(i)])
%     c = colorbar;
%     shading flat
% end
% 
%     subplot(3,2,5)
%     surface(x,y,N(:,:,1)')
%     title('N')
%     c = colorbar;
%     shading flat
%     
%     subplot(3,2,6)
%     surface(x,y,DOC(:,:,1)')
%     title('DOC')
%     c = colorbar;
%     shading flat
% 
% 
% 
