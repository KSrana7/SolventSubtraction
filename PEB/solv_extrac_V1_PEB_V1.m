% script for solvent subtraction: applied for the PEB data
% it contains:
%   Chemical Rank determination
%   Solvent Subtraction using Direct and orthogonal appraoches 


clc; clear all; close all;

set(groot,'defaultLineLineWidth',2)

addpath('K:\SolventSubtraction\ToolBoxes\tensor_toolbox-master\tensor_toolbox-master');
addpath('K:\SolventSubtraction\ToolBoxes\L-BFGS-B-C-master\L-BFGS-B-C-master\Matlab');
addpath('K:\SolventSubtraction\ToolBoxes\poblano_toolbox_1.1\poblano_toolbox_1.1');

addpath('K:\SolventSubtraction\ToolBoxes\MBtoolbox_v_02')
addpath('K:\SolventSubtraction\ToolBoxes\nway331')


%% Load Pre processed data

load 'Result/workspace_Data_Pre_Process_V0.mat'


%% Data arrangement and Determine Rank
global Z_ftir; global P_ftir;

Z_ftir=NaN(3,9,1764); 
% Data_ftir = combined_y;

Z_ftir(1,1:3,:)=Data_ftir(:,1:3)';
Z_ftir(2,4:6,:)=Data_ftir(:,4:6)';
Z_ftir(3,7:9,:)=Data_ftir(:,7:9)';



Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);


%% Use corcondia to det rank of FTIR
LOF_F=[];CF_F=[];
for i=1:8
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
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/Rank.png','-dpng','-r300')

% [A_f,B_f,C_f]=fac2let(Factors_ftir);
% M_f = nmodel(Factors_ftir);
%% 
R_f=3;


isequal(Z_ftir.*P_ftir, Z_ftir)


%% Run solvent substraction function
global N_I
N_I=50;

% tag='LG';
D_ftir_solvent=D_ftir_water_sampled;


dir_save = strcat('Result','/');

opts    = struct( 'factr', 1e-5, ...
        'pgtol', 0, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',500);

% method ='orthogonal';
method = 'DirectSubstraction';

Data_name = strjoin({'DD_',method,'.mat'},''); % dummy name for single file

Fun_SolvSub_V1(Z_ftir, P_ftir,Data_name,lam_ftir,D_ftir_solvent,dir_save, R_f,method,opts);


%% Move data for the Structure learning

sourceFile = fullfile(strjoin({'Result',method,strrep(Data_name,'.mat','.csv')},'\'));

destinationDir = fullfile(strjoin({'Label_input'},'\'));

if ~exist(destinationDir, 'dir')
    mkdir(destinationDir);
end

[~, fileName, fileExt] = fileparts(sourceFile);
destinationFile = fullfile(destinationDir, [fileName, fileExt]);

movefile(sourceFile, destinationFile);

fprintf('File moved to: %s\n', destinationFile);


%% to move without solvent substraction file
method = 'WO_Solv_Subs';

sourceFile = fullfile(strjoin({'Result',method,'DD_orthogonal_WOSolvExt_PCs.csv'},'\'));

destinationDir = fullfile(strjoin({'Label_input'},'\'));

if ~exist(destinationDir, 'dir')
    mkdir(destinationDir);
end

[~, fileName, fileExt] = fileparts(sourceFile);
destinationFile = fullfile(destinationDir, [fileName, fileExt]);

movefile(sourceFile, destinationFile);

fprintf('File moved to: %s\n', destinationFile);



