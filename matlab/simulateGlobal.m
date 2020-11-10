% Global run of size based unicelluar plankton model
% Transport Matrix used: MITgcm 2.8 degrees
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and be put into the location '../TMs'
%
% ../Cpp/model.cpp must be compiled locally to create model.so
%
% Input:
%  p: parameter structure from parametersGlobal
%  sim: (optional) simulation to use for initial conditions
%
function sim = simulateGlobal(p, sim) % if no initial conditions available: simulateGlobal(p)

% Load paths for switching TM
loadPath = '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
loadPath1 =  '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_';

%% Load Initial January TM:
load('../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_01.mat');
load('../TMs/MITgcm/grid.mat');
load('../TMs/MITgcm/config_data.mat');
load('../TMs/MITgcm/Matrix5/Data/boxes.mat');

% Preparing for timestepping. Using 43200 sec timesteps (0.5 days)
Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);

% number of elements in each layer:
lyr_n = permute(sum(sum(bathy)),[3,1,2]);
lyr_ind = [1 ; cumsum(lyr_n(1:end-1))+1];
lyr_end = cumsum(lyr_n);
%
% Load library:
%
%if bParallel
if isempty(gcp('nocreate'))
    parpool('AttachedFiles',{'../Cpp/model.so','../Cpp/model.h'});
end

h = gcp('nocreate');
poolsize = h.NumWorkers; 
parfor i=1:poolsize
    if libisloaded("model")
        unloadlibrary("model")
    end
    loadlibrary('../Cpp/model.so','../Cpp/model.h')
end
%
% Set parameters:
%
parfor i=1:poolsize
    calllib('model','setParameters', ...
        p.n, p.m, p.rhoCN, p.epsilonL, p.epsilonF, ...
        p.ANm, p.ALm, p.AFm, ...
        p.Jmax, p.JFmaxm, p.Jresp, p.Jloss_passive_m, ...
        p.theta, p.mort, p.mort2, p.mortHTL*p.mortHTLm, ...
        p.remin, p.remin2, p.cLeakage)
end
%%
% Initialize run:
%
simtime = p.tEnd*2; %simulation time in half days

month = 0;
iSave = 0;
YEAR = 0;

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
ticks = 1+2*cumsum(mon);
%
% Initial conditions:
%

if (nargin==2)
    disp('Starting from previous simulation.');
    N = double(sim.N(:,end));
    DOC = double(sim.DOC(:,end));
    Bmat = double(squeeze(sim.B(:,:,end)));
else
    load('../Data/N_oceandata.mat')
    DOC = zeros(nx,ny,nz) + p.DOC0;
    B = zeros(nx,ny,nz,p.n)+0.1; %biomass
    for i = 1:p.n
        B(:,:,1:5,i) = B(:,:,1:5,i)+p.B0(i);
    end
    % Convert from grid form to matrix form
    N   = gridToMatrix(N, [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    DOC = gridToMatrix(DOC, [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    Bmat = zeros(nb,p.n);
    for i = 1:p.n
        Bmat(:,i) = gridToMatrix(B(:,:,:,i), [], ...
            '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
    end
end
%
% Load environmental variables
%choose which temperature to use.

load('../TMs/MITgcm/BiogeochemData/Theta_bc.mat')

%Biochemical temperature
Tmat = zeros(nb,12);
for i = 1:12
    Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], '../TMs/MITgcm/Matrix5/Data/boxes.mat', '../TMs/MITgcm/grid.mat');
end
%
% Matrices for saving the solution:
%
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
Nm = single(zeros(nb,nSave));
DOCm = single(zeros(nb,nSave));
Bmatm= single(zeros(nb,p.n,nSave));
tSave = [];
%%
% Run transport matrix simulation
%
elapsed_time = zeros(1,simtime);
tic
for i=1:simtime
    telapsed = tic;
    %
    % Test for time to change monthly transport matrix
    %
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
            %disp(strcat(loadPath, num2str(month), '.mat'))
        else
            load(strcat(loadPath1, num2str(month), '.mat'));
            %disp(strcat(loadPath1, num2str(month), '.mat'))
        end
        
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
        
        % Set monthly mean temperature
        T = Tmat(:,month);
    end
    %%
    % Run Euler time step for half a day:
    %
    L = p.L(:,i);
    dt = p.dt;
    parfor k = 1:nb
        calllib('model','simulateEuler',[N(k); DOC(k); Bmat(k,:)'] ,...
            L(k), T(k),0,0, dt, 0.5);
    end
    if any(isnan([N;DOC;Bmat(:)]))
        warning('NaNs after running current grid box');
        keyboard
    end
    % Remove small concentrations
    N(N<1E-6) = 0;
    DOC(DOC<1E-6) = 0;
    Bmat(Bmat<1E-6) = 0;
    
    TIMESTEP = TIMESTEP+1;
    %
    % Transport
    %
    N   = Aimp * ( Aexp * N);
    DOC = Aimp * ( Aexp  * DOC);
    
    for k = 1:size(Bmat,2)
        Bmat(:,k) =  Aimp * (Aexp * Bmat(:,k));
    end
    
    elapsed_time(i) = toc(telapsed);
    %%
    % Save timeseries
    %
    if ((mod(i/2,p.tSave) < mod((i-1)/2,p.tSave)) || (i==simtime))
        iSave = iSave + 1;
        Nm(:,iSave) = (Nm(:,iSave)*(TIMESTEP-1) + single(N))/TIMESTEP;
        DOCm(:,iSave) = (DOCm(:,iSave)*(TIMESTEP-1) + single(DOC))/TIMESTEP;
        Bmatm(:,:,iSave)= (Bmatm(:,:,iSave)*(TIMESTEP-1) + single(Bmat))/TIMESTEP;
        tSave = [tSave, i*0.5];
        
        fprintf('t = %u days.\n',floor(i/2))
    end
end
time = toc;
fprintf('Solving time: %2u:%02u:%02u\n', ...
    [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
%
% Put results into sim structure:
%
sim.N = Nm;
sim.DOC = DOCm;
sim.B = Bmatm;
sim.t = tSave; % days where solution was saved
sim.p = p;
