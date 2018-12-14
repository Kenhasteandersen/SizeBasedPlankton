clear all
clc

% global C_L C_N C_F alpha_L alpha_N alpha_F phi_L phi_N phi_F M_L M_N M_F r_0 r_L r_N r_F beta_L beta_N beta_F C_CN m0 X_L X_F X_N Vlow1 Vhigh1 Vdiv1

%%%% Constant Environmental conditions %%%%
X_L1 = 250; % mu E/m^2/s From Maranon 2014 (Maranon et al 2013 Ecol Lett)
X_L = X_L1/4.6; % W/m^2
X_N = 70; 
X_F = 80; 

%%%% Constant investments %%%%
PhiL = 0.5;
PhiN = 0.2;
PhiF = 0.2;

%%%% Size limits %%%%
xlim_low=10^(-8);
xlim_high=10^(4);




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Affinities %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%_______________________
%%%%________ L ____________
%%%%_______________________

% subplot('Position',[ll 0.75 wd ht]) % subplot('Position',[left bottom width height])

nsamples=100000; % number of samples for Monte Carlo Simulation



%%%% data Taguchi %%%%

data_Taguchi = importdata('Taguchi.dat',',');
C=data_Taguchi(:,1); % pg C/cell
vol=data_Taguchi(:,2); % mu m^3
alpha=data_Taguchi(:,3); % mgC/(mg chlA) /h /(W m^-2)
CperChl=data_Taguchi(:,4); % mg C per chlA

xL04 = C; % pg/cell
xL4=10^(-6)*xL04; % converted to micro gc

% xL4= 1e-4*(3*vol/(4*pi)).^(1/3); % cm

C1 = C*10^(-6); % mug C/cell
chl = C1./(CperChl*10^3); % mug ChlA
alpha1 = alpha*10^(-3); % mug C/(mg chlA) /h /(W m^-2)
beta1 = alpha.*(chl); % mug C/h/(W m^-2)
beta2 = beta1*10^3; % mug C/h/(W m^-2)
beta3 = beta2*24; % mug C/day/(W m^-2)
yL4 = beta3; % mug C/day/(W m^-2)

loglog(xL4,yL4,'o','MarkerSize',5,'MarkerEdgeColor','k') % ,'MarkerFaceColor',[0,0,0]
hold on;



%%%% data Schwaderer 

data_Schwaderer = importdata('Edwards_2015_L_affinity_Data.csv',',');
% data_Schwaderer = importdata('Schwaderer_2011_light_affinity_data.csv',',');
xL0=data_Schwaderer(:,1); % given in micro m^3

xL01=10.^(0.74*log10(xL0)-0.58); % converted to pgc using the formula in Taguchi 1976
% xL01=10.^(0.76*log10(xL0)-0.29); % converted to pgc using the formula in Mullin et al 1966
% xL01=0.76*(xL0).^0.819; % (dinoflagellate) converted to pgc using the formula in Menden-Deuer and Lessard 2000
% xL01=0.228*(xL0).^0.811; % (diatom) converted to pgc using the formula in Menden-Deuer and Lessard 2000
xL=10^(-6)*xL01; % converted to mug C

yL1 = data_Schwaderer(:,2); % given in (mumol quanta)^-1 m^2 s day^-1
yL2 = 4.6*yL1; % W^-1 m^2 day^-1 (using 1W7m^2=4.6 mumol quanta/m^2-s)
yL = xL.*yL2; % mug C/day/(W m^-2)

loglog(xL,yL,'o','MarkerSize',5,'MarkerEdgeColor','k') % ,'MarkerFaceColor',[0,0,0]
hold on;



%%%%%%%% Our curve %%%%%%%%

alphaL = 1; %  10; % 0.5; %
CL = 0.02; %  0.01; % 0.05; %

AL=[];
AL1=[];
AL2=[];

V=[];
V1=[];

%%% First half of size %%%
xLlow=-8; 
xLhigh=-2; 
xLdiv=10^6;
xL1=logspace(xLlow,xLhigh,xLdiv);
V=xL1;
V1=V.^(2/3);

AL1 = (CL * V1.* V * PhiL * alphaL)./(PhiL * alphaL * V + CL * V1) ;

%%% Second half of size %%%
V=[];
V1=[];

xLlow=-2;
xLhigh=4;
xLdiv=10^6;
xL2=logspace(xLlow,xLhigh,xLdiv);
V=xL2;
V1=V.^(2/3);

AL2 = (CL * V1.* V * PhiL * alphaL)./(PhiL * alphaL * V + CL * V1) ;

%%% Combining both half %%%
xL(1:length(xL1))=xL1(:);
xL(length(xL1)+1:length(xL1)+length(xL2))=xL2(:);
V=xL;
AL(1:length(AL1))=AL1(:);
AL(length(AL1)+1:length(AL1)+length(AL2))=AL2(:);

loglog(V,AL,'k','LineWidth',2)

xlim([xlim_low xlim_high])
ylim([10^(-8) 10^(0)])
% xlabel('Size (\mug C)')
ylabel('\mug C/(W m^{-2} day)')
t = text(2*10^(-8),0.2,'(a) Affinity for light');
s = t.FontSize;
t.FontSize = 11;

ax = gca;
ax.XTick = [10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4)];
ax.YTick = [10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4)];
% set(gca,'xticklabel',{[]}) 
xlabel('\mug C')

