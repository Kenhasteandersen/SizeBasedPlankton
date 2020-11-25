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
% Output:
%  sim: structure with simulation results
%
function sim = simulateGlobal(p, sim) % if no initial conditions available: simulateGlobal(p)

% Load paths for switching TM
% loadPath = '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_0';
% loadPath1 =  '../TMs/MITgcm/Matrix5/TMs/matrix_nocorrection_';

%% Load Initial January TM:
load(strcat(p.pathMatrix0, '1.mat'));
load(p.pathGrid);
load(p.pathConfigData);
load(p.pathBoxes);

% Preparing for timestepping. Using 43200 sec timesteps (0.5 days)
Ix = speye(nb,nb);
Aexp = Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);
%
% Load library:
%
if p.bParallel
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',{'../Cpp/model.so','../Cpp/model.h'});
    end
   
    %
    % Set parameters:
    %
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        loadModel;
        setParameters(p);
    end
else
    loadModel;
    setParameters(p);
end

%
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
    load(p.pathN0)
    DOC = zeros(nx,ny,nz) + p.DOC0;
    B = zeros(nx,ny,nz,p.n); %biomass
    for i = 1:p.n
        B(:,:,1:2,i) = B(:,:,1:2,i)+p.B0(i);
    end
    % Convert from grid form to matrix form
    N   = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
    DOC = gridToMatrix(DOC, [], p.pathBoxes, p.pathGrid);
    Bmat = zeros(nb,p.n);
    for i = 1:p.n
        Bmat(:,i) = gridToMatrix(B(:,:,:,i), [], p.pathBoxes, p.pathGrid);
    end
end
%
% Load temperature:
%
load(p.pathTemp);
Tmat = zeros(nb,12);
for i = 1:12
    Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], p.pathBoxes, p.pathGrid);
end

%
% Load Light:
%

L0 = zeros(nb,730);
for i = 1:730
    L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-p.kw*Zbox);
end
%L(p.L<1) = 0;

%
% Matrices for saving the solution:
%
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
sim.N = single(zeros(nb,nSave));
sim.DOC = single(zeros(nb,nSave));
sim.B = single(zeros(nb,p.n,nSave));
sim.L = single(zeros(nb,nSave));
sim.T = single(zeros(nb,nSave));
dudt = zeros(1,2+p.n);
tSave = [];
%
% Run transport matrix simulation
%
elapsed_time = zeros(1,simtime);
disp('Starting simulation')
tic
for i=1:simtime
    telapsed = tic;
    %
    % Test for time to change monthly transport matrix
    %
    if ismember(i, ticks+730*YEAR)%
        month = month + 1;
        % Reset to Jan
        if month > 12
            month = 1;
        end
        % Load TM
        if month < 10
            load(strcat(p.pathMatrix0, num2str(month), '.mat'));
            %disp(strcat(p.pathMatrix0, num2str(month), '.mat'))
        else
            load(strcat(p.pathMatrix1, num2str(month), '.mat'));
            %disp(strcat(p.pathMatrix1, num2str(month), '.mat'))
        end
        
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
        
        % Set monthly mean temperature
        T = Tmat(:,month);
    end
    %
    % Run Euler time step for half a day:
    %
    L = L0(:,mod(i,365*2)+1);
    dt = p.dt;
    if p.bParallel
        parfor k = 1:nb
            u = [N(k); DOC(k); Bmat(k,:)'];
            u = calllib('model','simulateEuler', u, dudt, L(k), T(k),dt, 0.5);
            N(k) = max(0,u(1));
            DOC(k) = max(0,u(2));
            Bmat(k,:) = max(0,u(3:end))';
        end
    else
        for k = 1:nb
            u = [N(k); DOC(k); Bmat(k,:)'];
            u = calllib('model','simulateEuler', u, dudt, L(k), T(k), dt, 0.5);
            N(k) = u(1);
            DOC(k) = u(2);
            Bmat(k,:) = u(3:end)';
        end
    end
    if any(isnan([N;DOC;Bmat(:)]))
        warning('NaNs after running current grid box');
        keyboard
    end
    %
    % Transport
    %
    if p.bTransport
        N   = Aimp * ( Aexp * N);
        DOC = Aimp * ( Aexp  * DOC);
        
        for k = 1:size(Bmat,2)
            Bmat(:,k) =  Aimp * (Aexp * Bmat(:,k));
        end
    end
    elapsed_time(i) = toc(telapsed);
    %%
    % Save timeseries
    %
    if ((mod(i/2,p.tSave) < mod((i-1)/2,p.tSave)) || (i==simtime))
        iSave = iSave + 1;
        sim.N(:,iSave) = single(N);
        sim.DOC(:,iSave) = single(DOC);
        sim.B(:,:,iSave)= single(Bmat);
        sim.L(:,iSave) = L;
        sim.T(:,iSave) = T;
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
sim.t = tSave; % days where solution was saved
sim.p = p;

