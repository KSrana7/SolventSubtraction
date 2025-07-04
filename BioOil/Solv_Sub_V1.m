% Code of Solvent Subtraction: applied for the Bio-Oil data
%   Chemical Rank determination
%   Solvent Subtraction using Direct and orthogonal appraoches 

clc; clear all; close all;

set(groot,'defaultLineLineWidth',2)

addpath('K:\SolventSubtraction\ToolBoxes\tensor_toolbox-master\tensor_toolbox-master');
addpath('K:\SolventSubtraction\ToolBoxes\L-BFGS-B-C-master\L-BFGS-B-C-master\Matlab');
addpath('K:\SolventSubtraction\ToolBoxes\poblano_toolbox_1.1\poblano_toolbox_1.1');

addpath('K:\SolventSubtraction\ToolBoxes\MBtoolbox_v_02')
addpath('K:\SolventSubtraction\ToolBoxes\nway331')

%% Inputs


tag= 'WS';

files = dir(fullfile(pwd,strcat('\Result\preprocessed\',tag,'\','*.mat')));

% Data_name = 'Preprocessed_WIS 0_MeOH.mat';
files.name

Data_name = files(3).name%'Preprocessed_WIS 20_MeOH.mat';
solvent='water';% when O MeOH file is selected
solvent='20MeOH';%  when 1O and 20 MeOH file is selected




%% Load Data

files = dir(fullfile(pwd,strcat('\Result\preprocessed\',tag,'\','*.mat')));
index = find(strcmp({files.name}, Data_name));

filename_ftir = fullfile(pwd,strcat('\Result\preprocessed\',tag,'\',files(index).name));

loadedData = load(filename_ftir);

D_ftir = loadedData.Data_ftir_org;

lam_ftir=D_ftir(:,1);%wavenumber
Data_ftir=D_ftir(:,2:end); %intensity


% solvent data
filename_ftir_water = strcat('\Result\preprocessed\D_solv_water.mat');%'Aging reaction analysis/FTIR data/Pure water.CSV';
filename_ftir_MeOH = strcat('\Result\preprocessed\D_solv_MeOH.mat');%'Aging reaction analysis/FTIR data/MeOH.CSV';
filename_ftir_10_MeOH = strcat('\Result\preprocessed\D_solv_10MeOH.mat');%'Aging reaction analysis/FTIR data/10_ MeOH in water.CSV';
filename_ftir_20_MeOH = strcat('\Result\preprocessed\D_solv_20MeOH.mat');%'Aging reaction analysis/FTIR data/20_ MeOH in water.CSV';

loadedData = load(fullfile(pwd,filename_ftir_water));
Data_ftir_water= loadedData.Data_ftir_org_water;%xlsread(filename_ftir_water);  %data matrix
lam_ftir_water = Data_ftir_water(i_noise:end,1);
D_ftir_water = Data_ftir_water(i_noise:end,2);

loadedData = load(fullfile(pwd,filename_ftir_MeOH));
Data_ftir_MeOH= loadedData.Data_ftir_org_MeOH;%xlsread(filename_ftir_MeOH);  %data matrix
lam_ftir_MeOH = Data_ftir_MeOH(i_noise:end,1);
D_ftir_MeOH = Data_ftir_MeOH(i_noise:end,2);

loadedData = load(fullfile(pwd,filename_ftir_10_MeOH));
Data_ftir_10MeOH= loadedData.Data_ftir_org_10MeOH;%xlsread(filename_ftir_10_MeOH);  %data matrix
lam_ftir_10MeOH = Data_ftir_10MeOH(i_noise:end,1);
D_ftir_10MeOH = Data_ftir_10MeOH(i_noise:end,2);

loadedData = load(fullfile(pwd,filename_ftir_20_MeOH));
Data_ftir_20MeOH= loadedData.Data_ftir_org_20MeOH;%xlsread(filename_ftir_20_MeOH);  %data matrix
lam_ftir_20MeOH = Data_ftir_20MeOH(i_noise:end,1);
D_ftir_20MeOH = Data_ftir_20MeOH(i_noise:end,2);

if strcmp(solvent,'water')
    D_ftir_solvent=D_ftir_water;
elseif strcmp(solvent,'Methanol') || strcmp(solvent,'MeOH')
    D_ftir_solvent=D_ftir_MeOH;
elseif strcmp(solvent,'10Methanol') || strcmp(solvent,'10MeOH')
    D_ftir_solvent=D_ftir_10MeOH;
elseif strcmp(solvent,'20Methanol') || strcmp(solvent,'20MeOH')
    D_ftir_solvent=D_ftir_20MeOH;
else
    disp('Wrong solvent tag: Check again')
    return;
    
end

%% Determine Rank
global Z_ftir; global P_ftir;

Z_ftir=NaN(3,9,size(Data_ftir,1)); % time, Temp, wavenumber

Z_ftir(1,1:3,:)=Data_ftir(:,1:3)';
Z_ftir(2,4:6,:)=Data_ftir(:,4:6)';
Z_ftir(3,7:9,:)=Data_ftir(:,7:9)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);


% Use corcondia to det rank of FTIR
LOF_F=[];CF_F=[];
N_R=input('Enter Number of Parafac to run ');

for i=1:N_R
    [Factors_ftir,it_f,lof_f,cf_f]=parafac(Z_f,i);
    %M_f = nmodel(Factors_ftir);
    LOF_F=[LOF_F;lof_f];
    %cf_f = corcond(Z_f,Factors_ftir);
    CF_F=[CF_F;cf_f];
end

figure()
subplot(1,2,1)
plot(1:length(LOF_F),LOF_F,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Lack of fit (LOF)','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')
grid on
subplot(1,2,2)
plot(1:length(CF_F),CF_F,'-BX')
axis tight
xlabel('Number of components','fontweight','bold','FontSize',20)
ylabel('Core consistency','fontweight','bold','FontSize',20)
set(gca,'FontSize',20,'fontweight','bold')
grid on

fig = gcf;%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (13x7 inches, for example)
print(strcat('Result/plot/',tag,'Rank_ftir_',Data_name,'.png'),'-dpng','-r300')

% [A_f,B_f,C_f]=fac2let(Factors_ftir);
% M_f = nmodel(Factors_ftir);

% R_f=3;
R_f = input('Enter a number: ');

disp(['You entered: ', num2str(R_f)]);

isequal(Z_ftir.*P_ftir, Z_ftir)


%% Run solvent substraction function

opts    = struct( 'factr', 1e-5, ...
        'pgtol', 0, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',500);
dir_save = strcat('Result/',tag,'/');

% method ='orthogonal';
method = 'direct';
Fun_SolvSub_V1(Z_ftir, P_ftir,Data_name,lam_ftir,D_ftir_solvent,dir_save, R_f,method,opts);




