% script for solvent extraction: On Bitumen Data
% it contains:
%   Data generation that contains fraction of water spectra
%   Data preprocessing--> baseline correction and smoothening
%   Chemical Rank determination
%   Solvent Subtraction using Direct and orthogonal appraoches 
%   

clc; clear all; close all;

set(groot,'defaultLineLineWidth',2)


addpath('K:\SolventSubtraction\ToolBoxes\tensor_toolbox-master\tensor_toolbox-master');
addpath('K:\SolventSubtraction\ToolBoxes\L-BFGS-B-C-master\L-BFGS-B-C-master\Matlab');
addpath('K:\SolventSubtraction\ToolBoxes\poblano_toolbox_1.1\poblano_toolbox_1.1');

addpath('K:\SolventSubtraction\ToolBoxes\MBtoolbox_v_02')
addpath('K:\SolventSubtraction\ToolBoxes\nway331')



filename_ftir = 'Data/bitumenfinaldata.xlsx';
D_ftir= xlsread(filename_ftir);  %data matrix
T_ftir=D_ftir(1,2:end); %temperature
time_ftir=D_ftir(2,2:end); %time
X_ftir=D_ftir(3:end,:);%intensity
lam_ftir=X_ftir(:,1);%wavenumber
X11_ftir=X_ftir(:,2:end); %intensity
D_ftir=2-log10(X11_ftir);% change to absorbance which is interpreted as conc


filename_ftir_water = 'Data/FTIR_Water.xlsx';
Data_ftir_water= xlsread(filename_ftir_water);  %data matrix
lam_ftir_water = Data_ftir_water(:,1);
D_ftir_water = -log10(Data_ftir_water(:,2));

% subplot(2,1,1)
% plot(lam_ftir,D_ftir(:,1))
% subplot(2,1,2)
% plot(lam_ftir_water,0.05*D_ftir_water)


%% Data preprocessing
D_ftir_c = D_ftir;

%baseline and background correction
Dback_ftir=msbackadj(lam_ftir,D_ftir_c); % Correct baseline of signal--> spline approximation--> noise removal
Data_ftir=mssgolay(lam_ftir,Dback_ftir,'DEGREE',2,'Span',5);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal
Noise_ftir=Dback_ftir-Data_ftir;


% Plotting
%FTIR%
figure(1)
loc = round(linspace(1,42,6));
% labels = {'original','baseline adjusted','smoothening'};
for i=1:length(loc)
    subplot(2,3,i)
    plot(lam_ftir,D_ftir_c(:,loc(i)),'DisplayName','original'); hold on
    plot(lam_ftir,Dback_ftir(:,loc(i)),'DisplayName','baseline adjusted'); hold on
    plot(lam_ftir,Data_ftir(:,loc(i)),'DisplayName','smoothening')
    
    % Add legend
    legend('Location', 'best');
    
    % Add title
    title(sprintf("FTIR data Preprocessing:# %d",loc(i)));
    
    % Add axis labels
    xlabel('wavenumber');
    ylabel('Absorbance');
    grid on;
end

% Save figure 
% print('Result/plot/preprocessed_figure.png', '-dpng', '-r300');


%% merge signals

% Tolerance for accepting nearest x-values
tolerance = 0.5;

% Initialize combined signal
combined_lam = [];
combined_data = D_ftir;
D_ftir_water_sampled=zeros(length(lam_ftir),1);
total = [];

% Find the nearest x-values in signal 2 for each x-value in signal 1
for i = 1:length(lam_ftir)
    % Find the index of the nearest x-value
    [~, idx] = min(abs(lam_ftir_water - lam_ftir(i)));
    
    % Check if the nearest x-value is within tolerance
    if (abs(lam_ftir_water(idx) - lam_ftir(i))) <= tolerance
        combined_lam = [combined_lam, lam_ftir(i)];
        D_ftir_water_sampled(i) = 0.05*D_ftir_water(idx);
        combined_data(i,:) = combined_data(i,:) + D_ftir_water_sampled(i); % Add the corresponding y-values as per their weight fraction
    else
        % Print exception message
        total = [total, i];
        disp(['No matching x-value within tolerance for x = ', num2str(lam_ftir(i))]);
    end
end
%aaa = combined_data - D_ftir %validation
% Add Gaussian noise to the combined signal
noise_std = 0.0*mean(std(combined_data));%0.1; % Standard deviation for Gaussian noise% No noise in the original data
combined_y = combined_data + noise_std * randn(size(combined_data));

norm(combined_y-combined_data)

%% Data collection
D_ftir; % original ftir data
D_ftir_water_sampled; % water spectra downsampled to match the wavenumbeer
combined_data;% data with water spectra
combined_y; % data with water spectra and random noise 

% Plot the original and combined signals
figure;
loc = round(linspace(1,42,6));
for i=1:length(loc)
    subplot(2,3,i)
    plot(lam_ftir, D_ftir_water_sampled, 'r--',lam_ftir, D_ftir(:,loc(i)), 'b-',lam_ftir, combined_data(:,loc(i)));
    legend('water spectra', 'mixture spectra', 'Combined spectra', 'Location', 'best');
    xlabel('waveumber cm-1');
    ylabel('absorbance');
    title(sprintf('Combining Signal at (%d min, %d 0C)', [time_ftir(loc(i)),T_ftir(loc(i))]));
    grid on
end
% Save figure as PDF
% print('Result/plot/combined_figure.png', '-dpng', '-r300');



%% Determine Rank
global Z_ftir; global P_ftir;
Z_ftir=NaN(8,42,1764); 
Data_ftir = combined_y;

%Z_ftir(1,1,:)=Data_ftir(:,1)';
Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
Z_ftir(4,16,:)=Data_ftir(:,16)';
Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);


% Use corcondia to det rank of FTIR
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

fig = figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (13x7 inches, for example)
print('Result/plot/Rank_ftir.png','-dpng','-r300')

%--------------------------------------


R_f=4;  % based on the above plots


isequal(Z_ftir.*P_ftir, Z_ftir)


%% ===Data deconvolution on solvent free data=====================

Z_ftir=NaN(8,42,1764); 
Data_ftir = D_ftir;

%Z_ftir(1,1,:)=Data_ftir(:,1)';
Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
Z_ftir(4,16,:)=Data_ftir(:,16)';
Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);


% Initial guess
%[Ainit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,1)*tenmat(Z_ftir,1)'),R_f,[]);
%[Binit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,2)*tenmat(Z_ftir,2)'),R_f,[]);
%Cinit_ftir=double(tenmat(Z_ftir,3)*khatrirao(Binit_ftir,Ainit_ftir))*pinv((Binit_ftir'*Binit_ftir).*(Ainit_ftir'*Ainit_ftir));
%Minit_ftir={Ainit_ftir;Binit_ftir;Cinit_ftir};
%Minit_ftir = create_guess('Data',Z_ftir, 'Num_Factors', R_f, ...
  %    'Factor_Generator', 'nvecs');
%Call optimizer

SSE_ftir=[];M_FTIR = struct('tensrs',[]);

% opts    = struct( 'factr', 1e2, ...
%     'pgtol', 1e-9, ...
%     'm', 20, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',100);


for i=1:50
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%'opt_options',opts);%,'init', Minit_ftir);
SSE_ftir=[SSE_ftir;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

% csvwrite('Result/SSE_ftir.csv',SSE_ftir);

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir=M_FTIR(sse_locf).tensrs;


% csvwrite('Result/M_ftirU1.csv',M_ftir.U{1});
% csvwrite('Result/M_ftirU2.csv',M_ftir.U{2});
% csvwrite('Result/M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir.lambda, repmat(R1,1,ndims(M_ftir))); %<-- Create a diagonal core.
Y_ftir = ttensor(core1, M_ftir.U) ;%<-- Assemble the ttensor.


% error analysis
rmse_dd = sqrt(min(SSE_ftir));
T_A = reshape(M_ftir.u{3},[1764,R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_dd = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_dd)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_dd)])

% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir.u{2}(:,i));
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir.u{3}(:,i));
    legend('spectra','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution for solvent free data for factor:# %d",i));
    
end


figure()
plot(lam_ftir,Y_ftir.u{3});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none');
title(sprintf("FTIR data deconvolution for solvent free data: wavenumber vs Spectra"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('Baseline Values \n RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_dd,Mcos_l_dd]), ...
    'FitBoxToText', 'on');


% Saving
current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_solventFreedata_R_4_', current_timestamp, '.mat'];
save(file_name, 'M_ftir','Y_ftir','lam_ftir','D_ftir_water_sampled');

% 
fig = figure(2);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Solv_free/DD_solv_free_R_4.png','-dpng','-r300')

%% ===Data deconvolution Without Solvent Subtraction=================================
% Fit a PARAFAC model for Z_ftir

Z_ftir=NaN(8,42,1764); 
Data_ftir = combined_y;

%Z_ftir(1,1,:)=Data_ftir(:,1)';
Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
Z_ftir(4,16,:)=Data_ftir(:,16)';
Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);

SSE_ftir_ws=[];M_FTIR_ws = struct('tensrs',[]);

for i=1:50  % best of 50 repetitions 
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%,'init', Minit_ftir);
SSE_ftir_ws=[SSE_ftir_ws;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR_ws(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

[sse_minf,sse_locf]=min(SSE_ftir_ws)

M_ftir_ws=M_FTIR_ws(sse_locf).tensrs;

% 
% csvwrite('Result/M_ftirU1.csv',M_ftir.U{1});
% csvwrite('Result/M_ftirU2.csv',M_ftir.U{2});
% csvwrite('Result/M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir_ws.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_ws.lambda, repmat(R1,1,ndims(M_ftir_ws))); %<-- Create a diagonal core.
Y_ftir_ws = ttensor(core1, M_ftir_ws.U) ;%<-- Assemble the ttensor.



% error analysis
rmse_ws = sqrt(min(SSE_ftir_ws));
T_A = reshape(M_ftir_ws.u{3},[1764,R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])

% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir_ws.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir_ws.u{2}(:,i));
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_ws.u{3}(:,i));
    legend('spectra','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
    
end

figure()
org_loc = [1,2,3,4];
for i=1:R1
    
    subplot(R1,1,i)
    plot(lam_ftir,[Y_ftir_ws.u{3}(:,i),Y_ftir.u{3}(:,org_loc(i))]);
    legend(sprintf('PC_%s',num2str(i)),'Solv_Free_PC','interpreter','none','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
    
end

fig = figure(4);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/WO_Solv_Ext/DD_wo_solv_ext_PCs_R_4.png','-dpng','-r300')


figure()
plot(lam_ftir,Y_ftir_ws.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
title(sprintf("FTIR data deconvolution without solvent extraction: wavenumber vs Spectra"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_ws,Mcos_l_ws]), ...
    'FitBoxToText', 'on');



%Saving
fig = figure(5);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/WO_Solv_Ext/DD_wo_solv_ext_R_4.png','-dpng','-r300')



current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_WO_solEx_R_4_', current_timestamp, '.mat'];
save(file_name, 'M_ftir_ws','Y_ftir_ws','lam_ftir','D_ftir_water_sampled');






P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);


% Initial guess
%[Ainit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,1)*tenmat(Z_ftir,1)'),R_f,[]);
%[Binit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,2)*tenmat(Z_ftir,2)'),R_f,[]);
%Cinit_ftir=double(tenmat(Z_ftir,3)*khatrirao(Binit_ftir,Ainit_ftir))*pinv((Binit_ftir'*Binit_ftir).*(Ainit_ftir'*Ainit_ftir));
%Minit_ftir={Ainit_ftir;Binit_ftir;Cinit_ftir};
%Minit_ftir = create_guess('Data',Z_ftir, 'Num_Factors', R_f, ...
  %    'Factor_Generator', 'nvecs');
%Call optimizer

SSE_ftir=[];M_FTIR = struct('tensrs',[]);

% opts    = struct( 'factr', 1e2, ...
%     'pgtol', 1e-9, ...
%     'm', 20, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',100);


for i=1:50
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%'opt_options',opts);%,'init', Minit_ftir);
SSE_ftir=[SSE_ftir;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

% csvwrite('Result/SSE_ftir.csv',SSE_ftir);

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir=M_FTIR(sse_locf).tensrs;


% csvwrite('Result/M_ftirU1.csv',M_ftir.U{1});
% csvwrite('Result/M_ftirU2.csv',M_ftir.U{2});
% csvwrite('Result/M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir.lambda, repmat(R1,1,ndims(M_ftir))); %<-- Create a diagonal core.
Y_ftir = ttensor(core1, M_ftir.U) ;%<-- Assemble the ttensor.


% error analysis
rmse_dd = sqrt(min(SSE_ftir));
T_A = reshape(M_ftir.u{3},[1764,R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_dd = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_dd)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_dd)])

% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir.u{2}(:,i));
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir.u{3}(:,i));
    legend('spectra','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution for solvent free data for factor:# %d",i));
    
end


figure()
plot(lam_ftir,Y_ftir.u{3});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none');
title(sprintf("FTIR data deconvolution for solvent free data: wavenumber vs Spectra"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('Baseline Values \n RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_dd,Mcos_l_dd]), ...
    'FitBoxToText', 'on');


% Saving
current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_solventFreedata_R_4_', current_timestamp, '.mat'];
save(file_name, 'M_ftir','Y_ftir','lam_ftir','D_ftir_water_sampled');

% 
fig = figure(2);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Solv_free/DD_solv_free_R_4.png','-dpng','-r300')

%% ===Data deconvolution Without Solvent Subtraction=================================
% Fit a PARAFAC model for Z_ftir

Z_ftir=NaN(8,42,1764); 
Data_ftir = combined_y;

%Z_ftir(1,1,:)=Data_ftir(:,1)';
Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
Z_ftir(4,16,:)=Data_ftir(:,16)';
Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);

SSE_ftir_ws=[];M_FTIR_ws = struct('tensrs',[]);

for i=1:50  % best of 50 repetitions 
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%,'init', Minit_ftir);
SSE_ftir_ws=[SSE_ftir_ws;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR_ws(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

[sse_minf,sse_locf]=min(SSE_ftir_ws)

M_ftir_ws=M_FTIR_ws(sse_locf).tensrs;

% 
% csvwrite('Result/M_ftirU1.csv',M_ftir.U{1});
% csvwrite('Result/M_ftirU2.csv',M_ftir.U{2});
% csvwrite('Result/M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir_ws.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_ws.lambda, repmat(R1,1,ndims(M_ftir_ws))); %<-- Create a diagonal core.
Y_ftir_ws = ttensor(core1, M_ftir_ws.U) ;%<-- Assemble the ttensor.



% error analysis
rmse_ws = sqrt(min(SSE_ftir_ws));
T_A = reshape(M_ftir_ws.u{3},[1764,R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])

% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir_ws.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir_ws.u{2}(:,i));
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_ws.u{3}(:,i));
    legend('spectra','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
    
end

figure()
org_loc = [1,2,3,4];
for i=1:R1
    
    subplot(R1,1,i)
    plot(lam_ftir,[Y_ftir_ws.u{3}(:,i),Y_ftir.u{3}(:,org_loc(i))]);
    legend(sprintf('PC_%s',num2str(i)),'Solv_Free_PC','interpreter','none','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
    
end

fig = figure(4);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/WO_Solv_Ext/DD_wo_solv_ext_PCs_R_4.png','-dpng','-r300')


figure()
plot(lam_ftir,Y_ftir_ws.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
title(sprintf("FTIR data deconvolution without solvent extraction: wavenumber vs Spectra"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_ws,Mcos_l_ws]), ...
    'FitBoxToText', 'on');



%Saving
fig = figure(5);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/WO_Solv_Ext/DD_wo_solv_ext_R_4.png','-dpng','-r300')



current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_WO_solEx_R_4_', current_timestamp, '.mat'];
save(file_name, 'M_ftir_ws','Y_ftir_ws','lam_ftir','D_ftir_water_sampled');





%% ===Data deconvolution With Solvent Subtraction=================================

%  Direct approach

SSE_ftir_s=[];M_FTIR_s = struct('tensrs',[]);


%solvent information
%if solvent spectra is known
Ad = 1.0*D_ftir_water_sampled;%==================Signature spectra (Original)

S_O={};


%%_________________________________________________________________
% First only solvent substraction from the tensor data, later data
% deconvolution on solvent free tensor data
%_________________________________________________________________
% *****1.first solvent substraction and later latent data deconvolution
SSE_ftir=[];

for i=1:1
opt_options    = struct( 'factr', 1e0, ...
    'pgtol', 1e-9, ...
    'm', 10, ...
    'maxIts',10000, ...
    'maxTotalIts',50000, ...
    'printEvery',10);

R2 = 1;
[M_ftir1,~,S_O{i},output_ftir1] = cp_wopt_s22(Z_ftir, P_ftir, R2,Ad,'lower',0,'opt_options',opt_options);%,'init', M_ftir.u);

Z_ftir_s = Z_ftir -tensor(S_O{i});
% 
opts    = struct( 'factr', 1e0, ...
    'pgtol', 1e-9, ...
    'm', 10, ...
    'maxIts',10000, ...
    'maxTotalIts',50000, ...
    'printEvery',500);
    
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir_s, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  

% SSE_ftir_s=[1000;sum(sum(sum((double(Z_ftir-tensor(arrange(S_O{i}))-tensor(M_ftir))).^2)))];
% j=1;
% 
% while 100*abs((SSE_ftir_s(j)-SSE_ftir_s(j+1))/SSE_ftir_s(j))>0.1 
%     
%     Z_ftir_s1 = Z_ftir -tensor(M_ftir);
%     [M_ftir1,~,S_O{i},output_ftir] = cp_wopt_s22(Z_ftir_s1, P_ftir, R_f,Ad,'lower',0,'init',{S_O{i}{1};S_O{i}{2}});
%     
%     Z_ftir_s2 = Z_ftir -tensor(arrange(S_O{i}));
%     
%     [M_ftir,~,output_ftir] = cp_wopt(Z_ftir_s2, P_ftir, R_f,'lower',0,'init', M_ftir);  
%     SSE_ftir_s=[SSE_ftir_s;sum(sum(sum((double(Z_ftir-tensor(arrange(S_O{i}))-tensor(M_ftir))).^2)))];
% 
%     j=j+1;
% 
%     disp(['Iteration for ALS: ',num2str(j)])
% 
% end
SSE_ftir=[SSE_ftir;sum(sum(sum((double(tensor(M_ftir)-P_ftir.*Z_ftir_s)).^2)))];
M_FTIR_s(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end
% ****1
%_________________________________________________________________


[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir_s=M_FTIR_s(sse_locf).tensrs;
% S_Out = S_O{sse_locf};


%convert to a ttensor
R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.


% error analysis

rmse_ws = sqrt(min(SSE_ftir));
T_A = reshape(M_ftir_s.u{3}(:,1:R_f),[size(M_ftir_s.u{3},1),R_f]);
%-----------------------------------------------------%
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])


% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir_s.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir_s.u{2}(:,i)); hold on
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_s.u{3}(:,i))
    legend('spectra','Interpreter', 'none','Location', 'best');
    sgtitle(sprintf("FTIR data deconvolution for factor:# %d",i));
    grid on;

end


% load 'Result/workspace_DD_solventFreedata_2024-02-26_12-30'
load 'Result/workspace_DD_WO_solEx_R_4_2025-06-30_15-39'

figure()
org_loc = [1,3,2,4];
for i=1:R_f
    subplot(R_f,1,i)
    plot(lam_ftir,[Y_ftir_s.u{3}(:,i),Y_ftir.u{3}(:,org_loc(i))]);
    legend(sprintf('PC_%s',num2str(i)),'Solv_free_PC','interpreter','none','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
  
end



%Saving
fig = figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Direct/DD_solv_ext_PCs_R_3.png','-dpng','-r300')



figure()
subplot(3,1,1)
plot(lam_ftir,M_ftir_init{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
title(sprintf("(Initialisation)"));
grid on;

subplot(3,1,2)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(3,1,3)
plot(lam_ftir,Y_ftir_s.u{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
grid on;


%Saving
fig = figure(2);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Direct/DD_solv_ext_All_R_3.png','-dpng','-r300')





figure()
subplot(2,1,1)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(2,1,2)
% Y2 = tensor(S_Out);%P_ftir.*full(ktensor(S_Out));
S_Out = {Y_ftir_s.u{1}(:,R_f+1);Y_ftir_s.u{2}(:,R_f+1);Y_ftir_s.u{3}(:,R_f+1)};
Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;

mean_tensor = zeros(size(Y2,3),1);
for i=1:size(Y2,3)
    Y_ext = double(Y2(:,:,i));
    Y_ext(Y_ext == 0) = NaN;
    mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
end
plot(1:length(Ad),[Ad,mean_tensor]);% when solvent profile is known
legend('Ground truth','Extracted')
sgtitle(sprintf("Solvent Substraction"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similaity:%.2f', [rmse_ws,Mcos_l_ws]), ...
    'FitBoxToText', 'on');



%Saving
fig = figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/DD_solv_ext_first_PCs.png','-dpng','-r300')

current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_solEx_', current_timestamp, '.mat'];
save(file_name, 'M_ftir_s','Y_ftir_s','S_Out','lam_ftir','D_ftir_water_sampled');

%% ===Data deconvolution With Solvent Subtraction=================================

%  Orthogonal approach

SSE_ftir_s=[];M_FTIR_s = struct('tensrs',[]);


%solvent information
%if solvent spectra is known
Ad = 1.0*D_ftir_water_sampled;%==================Signature spectra (Original)

S_O={};

% for i=1:1
% [M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%,'init', Minit_ftir);  
% % [M_ftir,~,S_O{i},output_ftir] = cp_wopt_s2(Z_ftir, P_ftir, R_f,Ad,'lower',0);%,'init', Minit_ftir);
% [M_ftir,~,S_O{i},output_ftir] = cp_wopt_s21(Z_ftir, P_ftir, R_f,Ad,'lower',0,'init', M_ftir.u);
% SSE_ftir_s=[SSE_ftir_s;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
% M_FTIR_s(i).tensrs=M_ftir;
% disp(['iteration, ' num2str(i)])
% end


%_________________________________________________________________
% Orthogonal CP decomposition
%_________________________________________________________________
SSE_ftir=[];

for i=1:50

opts    = struct( 'factr', 1e-3, ...
        'pgtol', 1e-5, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',500);
   
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  
% 
% [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
opts    = struct( 'factr', 1e-3, ...
    'pgtol', 0e-5, ...
    'm', 5, ...
    'maxIts',50000, ...
    'maxTotalIts',50000, ...
    'printEvery',100);

M_ftir_init = M_ftir.u;
for j=1:size(M_ftir_init,1)
    if j==3
        M_ftir_init{j} = horzcat(M_ftir_init{j},abs(Ad)); 
    else
        M_ftir_init{j} = horzcat(M_ftir_init{j},abs(matrandnorm(size(M_ftir_init{j},1),1))); 
        
    end
end

% Get solvent free spectra
[M_ftir,~,mu_l_o,output_ftir] = cp_wopt_s26(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts,'init', M_ftir_init);

% Fix spectra and recalibrate the concentration
A_s = M_ftir.u{3};
R_c = R_f+1;
[M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}},'verbosity',0);

SSE_ftir=[SSE_ftir;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];

M_FTIR_s(i).tensrs=M_ftir_c;
disp(['iteration, ' num2str(i)])
disp(output_ftir)
end
% 

% figure();
% plot([output_ftir.OptOut.err]);
% grid on;
% legend('func','grad')


%_________________________________________________________________

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir_s=M_FTIR_s(sse_locf).tensrs;
% S_Out = S_O{sse_locf};


%convert to a ttensor
R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.


% error analysis

M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s.u{1}(:,1:R_f);M_ftir_s.u{2}(:,1:R_f);M_ftir_s.u{3}(:,1:R_f)});
rmse_ws = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
T_A = reshape(M_ftir_s.u{3}(:,1:R_f),[size(M_ftir_s.u{3},1),R_f]);

profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])


% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir_s.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir_s.u{2}(:,i)); hold on
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_s.u{3}(:,i))
    legend('spectra','Interpreter', 'none','Location', 'best');
    sgtitle(sprintf("FTIR data deconvolution for factor:# %d",i));
    grid on;

end


% load 'Result/workspace_DD_solventFreedata_2024-02-26_12-30'
load 'Result/workspace_DD_WO_solEx_R_4_2025-06-30_15-39'

figure()
org_loc = [1,3,2,4];
for i=1:R_f
    subplot(R_f,1,i)
    plot(lam_ftir,[Y_ftir_s.u{3}(:,i),Y_ftir.u{3}(:,org_loc(i))]);
    legend(sprintf('PC_%s',num2str(i)),'Solv_free_PC','interpreter','none','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
  
end



%Saving
fig = figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Orthogonal/DD_solv_ext_PCs_R_3.png','-dpng','-r300')



figure()
subplot(3,1,1)
plot(lam_ftir,M_ftir_init{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
title(sprintf("(Initialisation)"));
grid on;

subplot(3,1,2)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(3,1,3)
plot(lam_ftir,Y_ftir_s.u{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
grid on;


%Saving
fig = figure(2);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Orthogonal/DD_solv_ext_All_R_3.png','-dpng','-r300')





figure()
subplot(2,1,1)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(2,1,2)
% Y2 = tensor(S_Out);%P_ftir.*full(ktensor(S_Out));
S_Out = {Y_ftir_s.u{1}(:,R_f+1);Y_ftir_s.u{2}(:,R_f+1);Y_ftir_s.u{3}(:,R_f+1)};
Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;

mean_tensor = zeros(size(Y2,3),1);
for i=1:size(Y2,3)
    Y_ext = double(Y2(:,:,i));
    Y_ext(Y_ext == 0) = NaN;
    mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
end
plot(1:length(Ad),[Ad,mean_tensor]);% when solvent profile is known
legend('Ground truth','Extracted')
sgtitle(sprintf("Solvent Substraction"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similaity:%.2f', [rmse_ws,Mcos_l_ws]), ...
    'FitBoxToText', 'on');



%Saving
fig = figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/DD_solv_ext_first_PCs.png','-dpng','-r300')

current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
file_name = ['Result/workspace_DD_solEx_', current_timestamp, '.mat'];
save(file_name, 'M_ftir_s','Y_ftir_s','S_Out','lam_ftir','D_ftir_water_sampled');



%% ===USE CASE I: solvent prior information : partial--> Solvent concentration is not exactly known=======

for j=[0.01,0.1,0.5,1.5,2.0,2.5,5.0]%[0.01,0.1,0.5,1.5,2.0,2.5,5.0]

SSE_ftir_s2=[];M_FTIR_s2 = struct('tensrs',[]);

% %solvent information
% %if solvent spectra is known
Ad = D_ftir_water_sampled;

%______________________First solvent substraction, later
%deconvolution________
% S_O_2={};A_O={};
% for i=1:50
%     opt_options    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',10);
% 
%     R2 =1;
% 
%     [M_ftir1,~,S_O_2{i},output_ftir1] = cp_wopt_s22(Z_ftir, P_ftir, R2,j*Ad,'lower',0,'opt_options',opt_options);%,'init', M_ftir.u);
% 
%     Z_ftir_s = Z_ftir -tensor(S_O_2{i});
% 
%     opts    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',500);
% 
%     [M_ftir,~,output_ftir] = cp_wopt(Z_ftir_s, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  
% 
%     % [M_ftir,~,S_O_2{i},output_ftir] = cp_wopt_s2(Z_ftir, P_ftir, R_f,j*D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
% %     [M_ftir,~,S_O_2{i},A_O{i},output_ftir] = cp_wopt_s3(Z_ftir, P_ftir, R_f,D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
%     SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(tensor(M_ftir)-P_ftir.*Z_ftir_s)).^2)))];
%     M_FTIR_s2(i).tensrs=M_ftir;
%     disp(['iteration, ' num2str(i)])
%     
% end


%_________________________________________________________________
% Orthogonal CP decomposition
%_________________________________________________________________
% SSE_ftir=[];
disp(['Calculation for ',num2str(j)])

for i=1:50

opts    = struct( 'factr', 1e-3, ...
        'pgtol', 1e-5, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',0);
   
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0,'opt_options',opts,'verbosity',0);%,'init', Minit_ftir);  
% 
% [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
opts    = struct( 'factr', 1e-3, ...
    'pgtol', 0e-5, ...
    'm', 5, ...
    'maxIts',50000, ...
    'maxTotalIts',50000, ...
    'printEvery',0);

M_ftir_init = M_ftir.u;
for jj=1:size(M_ftir_init,1)
    if jj==3
        M_ftir_init{jj} = horzcat(M_ftir_init{jj},abs(j*Ad)); 
    else
        M_ftir_init{jj} = horzcat(M_ftir_init{jj},abs(matrandnorm(size(M_ftir_init{jj},1),1))); 
        
    end
end

% Get solvent Free Spectra
[M_ftir,~,mu_l_o,output_ftir_s] = cp_wopt_s26(Z_ftir, P_ftir, R_f,j*Ad,'lower',0,'opt_options',opts,'init', M_ftir_init,'verbosity',0);

% Fix spectra and recalibrate the concentration
A_s = M_ftir.u{3};
R_c = R_f+1;
[M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}},'verbosity',0);

% figure();plot(output_ftir_c.OptOut.err); grid on;
% sum(sum(sum((double(P_ftir.*tensor(M_ftir)-Z_ftir)).^2)))

SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];
M_FTIR_s2(i).tensrs=M_ftir_c;
% disp(['iteration, ' num2str(i)])
% disp(output_ftir_c)
end
% 

disp('Data deconvolution for intial guess')
disp(output_ftir)
disp('Data deconvolution after solvent extraction')
disp(output_ftir_s)
disp('Data deconvolution after calibration')
disp(output_ftir_c)
    
% % % csvwrite('Result/SSE_ftir_s.csv',SSE_ftir_s);

[sse_minf,sse_locf]=min(SSE_ftir_s2);

M_ftir_s2=M_FTIR_s2(sse_locf).tensrs;
% S_Out2 = S_O_2{sse_locf};

% error analysis
M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s2.u{1}(:,1:R_f);M_ftir_s2.u{2}(:,1:R_f);M_ftir_s2.u{3}(:,1:R_f)});
rmse_s2 = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
T_A = reshape(M_ftir_s2.u{3}(:,1:R_f),[size(M_ftir_s2.u{3},1),R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_s2 = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
% disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_s2)])
% disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_s2)])


% csvwrite('Result/M_ftirU1_s.csv',M_ftir_s.U{1});
% csvwrite('Result/M_ftirU2_s.csv',M_ftir_s.U{2});
% csvwrite('Result/M_ftirU3_s.csv',M_ftir_s.U{3});

%convert to a ttensor
R1 = length(M_ftir_s2.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s2.lambda, repmat(R1,1,ndims(M_ftir_s2))); %<-- Create a diagonal core.
Y_ftir_s2 = ttensor(core1, M_ftir_s2.U) ;%<-- Assemble the ttensor.

% 
% R1 = length(S_Out2.lambda);  %<-- Number of factors in X.
% core1 = tendiag(S_Out2.lambda, repmat(R1,1,ndims(S_Out2))); %<-- Create a diagonal core.
% Y_ftir_s1 = ttensor(core1, S_Out2.U) ;%<-- Assemble the ttensor.
% 

figure()
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;

subplot(2,1,1)
plot(lam_ftir,Y_ftir_s2.u{3}(:,1:R_f),lam_ftir, j*Ad,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra for %.2f",j));
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_s2,Mcos_l_s2]), ...
    'FitBoxToText', 'on');
grid on;

subplot(2,1,2)
% Y2 = tensor(S_Out2);%P_ftir.*full(ktensor(S_Out));
Y2= P_ftir.*ktensor({M_ftir_s2.u{1}(:,R_f+1);M_ftir_s2.u{2}(:,R_f+1);M_ftir_s2.u{3}(:,R_f+1)});
% S_Out = {Y_ftir_s2.u{1}(:,R_f+1);Y_ftir_s2.u{2}(:,R_f+1);Y_ftir_s2.u{3}(:,R_f+1)};
% Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;
mean_tensor = zeros(size(Y2,3),1);
for i=1:size(Y2,3)
    Y_ext = double(Y2(:,:,i));
    Y_ext(Y_ext == 0) = NaN;
    mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
end
plot(1:length(D_ftir_water_sampled),[D_ftir_water_sampled,mean_tensor])
legend('Ground truth','Mean_Extracted','Interpreter', 'none','Location', 'best')
sgtitle(sprintf("Solvent Substraction for offset factor %.2f",j));
grid on;

% Adjust figure size
fig = gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)

print(sprintf('Result/plot/Orthogonal/Int_C_var/DD_SolvExt_orth_4_times_%s.png',num2str(j)),'-dpng','-r300')

end
% saving
% 
% figure(13)
% print('Result/plot/Solv_Ext_full/DD_SolvExt_full_times_1_4.png','-dpng','-r300')

% current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
% file_name = ['Result/workspace_DD_solEx_', current_timestamp, '.mat'];
% save(file_name, 'M_ftir_s2','Y_ftir_s2','S_Out2','lam_ftir','D_ftir_water_sampled');


%% === USE CASE II: Performance in presenace of noise============================

jj=[0.0,0.1,0.5,1.5,2.0];S_Out2={};SNR=[];
for j =1:length(jj)%[0.01,0.1,0.5,1.5,2.0]
% Add Gaussian noise to the combined signal
noise_std = jj(j)*mean(std(combined_data));%0.1; % Standard deviation for Gaussian noise
combined_y = combined_data + noise_std * randn(size(combined_data));

SNR(j) = mean(mean(combined_y)./std(combined_y));

sprintf('SNR for J %d, %d',SNR(j),j)


Z_ftir=NaN(8,42,1764); 
Data_ftir = combined_y;

%Z_ftir(1,1,:)=Data_ftir(:,1)';
Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
Z_ftir(4,16,:)=Data_ftir(:,16)';
Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
Z_f=Z_ftir;
Z_ftir(isnan(Z_ftir))=0;

P_ftir=(Z_ftir~=0);
P_ftir=tensor(P_ftir);
    
Z_ftir=tensor(Z_ftir);



SSE_ftir_s=[];SSE_ftir_s2=[];M_FTIR_s2 = struct('tensrs',[]);

% %solvent information
% %if solvent spectra is known
Ad = D_ftir_water_sampled;

S_O_2={};A_O={};

%_________________________first solvent extraction, later deconvolution____
% for i=1:50
%     opt_options    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',10);
% 
%     R2 = 2;
%     [M_ftir1,~,S_O_2{i},output_ftir] = cp_wopt_s22(Z_ftir, P_ftir, R2,Ad,'lower',0,'opt_options',opt_options);%,'init', M_ftir.u);
% 
%     Z_ftir_s = Z_ftir -tensor(S_O_2{i});
% 
%     opts    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',500);
% 
%     [M_ftir,~,output_ftir2] = cp_wopt(Z_ftir_s, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  
% 
%     % [M_ftir,~,S_O_2{i},output_ftir] = cp_wopt_s2(Z_ftir, P_ftir, R_f,j*D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
% %     [M_ftir,~,S_O_2{i},A_O{i},output_ftir] = cp_wopt_s3(Z_ftir, P_ftir, R_f,D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
%     SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(tensor(M_ftir)-P_ftir.*Z_ftir_s)).^2)))];
%     M_FTIR_s2(i).tensrs=M_ftir;
%     disp(['iteration, ' num2str(i)])
% end
%_____________________________________________________________________________

%________________________________________________________________________________________________
%___________________________Orthogonal deconvolution with solvent as additional factor________
% ____________________________________________________________________________________________


for i=1:50

opts    = struct( 'factr', 1e-3, ...
        'pgtol', 1e-5, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',500);
   
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  
% 
% [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
opts    = struct( 'factr', 1e-3, ...
    'pgtol', 0e-5, ...
    'm', 5, ...
    'maxIts',50000, ...
    'maxTotalIts',50000, ...
    'printEvery',100);

M_ftir_init = M_ftir.u;
for ij=1:size(M_ftir_init,1)
    if ij==3
        M_ftir_init{ij} = horzcat(M_ftir_init{ij},abs(j*Ad)); 
    else
        M_ftir_init{ij} = horzcat(M_ftir_init{ij},abs(matrandnorm(size(M_ftir_init{ij},1),1))); 
        
    end
end

% Get solvent Free Spctra
[M_ftir,~,mu_l_o,output_ftir] = cp_wopt_s26(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts,'init', M_ftir_init);

% Fix spectra and recalibrate the concentration
A_s = M_ftir.u{3};
R_c = R_f+1;
[M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}});

% figure();plot(output_ftir_c.OptOut.err); grid on;
% sum(sum(sum((double(P_ftir.*tensor(M_ftir)-Z_ftir)).^2)))

SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];
M_FTIR_s2(i).tensrs=M_ftir_c;
disp(['iteration, ' num2str(i)])
disp(output_ftir_c)
end
%________________________________________________________________

% % % csvwrite('Result/SSE_ftir_s.csv',SSE_ftir_s);

[sse_minf,sse_locf]=min(SSE_ftir_s2);

M_ftir_s2=M_FTIR_s2(sse_locf).tensrs;
% S_Out2{j} = S_O_2{sse_locf};

% error analysis
M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s2.u{1}(:,1:R_f);M_ftir_s2.u{2}(:,1:R_f);M_ftir_s2.u{3}(:,1:R_f)});
rmse_s2 = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
T_A = reshape(M_ftir_s2.u{3}(:,1:R_f),[size(M_ftir_s2.u{3},1),R_f]);
% rmse_s2 = sqrt(min(SSE_ftir_s2));
% T_A = reshape(M_ftir_s2.u{3}(:,1:R_f),[1764,4]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_s2 = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_s2)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_s2)])


%convert to a ttensor
R1 = length(M_ftir_s2.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s2.lambda, repmat(R1,1,ndims(M_ftir_s2))); %<-- Create a diagonal core.
Y_ftir_s2 = ttensor(core1, M_ftir_s2.U) ;%<-- Assemble the ttensor.


figure()
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;

subplot(2,1,1)
plot(lam_ftir,5*Y_ftir_s2.u{3}(:,1:R_f),lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra(SNR= %.2f)",SNR(j)));
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_s2,Mcos_l_s2]), ...
    'FitBoxToText', 'on');
grid on;

subplot(2,1,2)
% Y2 = tensor(S_Out2);%P_ftir.*full(ktensor(S_Out));
Y2= P_ftir.*ktensor({M_ftir_s2.u{1}(:,R_f+1);M_ftir_s2.u{2}(:,R_f+1);M_ftir_s2.u{3}(:,R_f+1)});
% S_Out = {Y_ftir_s2.u{1}(:,R_f+1);Y_ftir_s2.u{2}(:,R_f+1);Y_ftir_s2.u{3}(:,R_f+1)};
% Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;
mean_tensor = zeros(size(Y2,3),1);
for i=1:size(Y2,3)
    Y_ext = double(Y2(:,:,i));
    Y_ext(Y_ext == 0) = NaN;
    mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
end
plot(1:length(D_ftir_water_sampled),[D_ftir_water_sampled,mean_tensor])
legend('Ground truth','Mean_Extracted','Interpreter', 'none','Location', 'best')
title(sprintf("Solvent Substraction,SNR= %.2f",SNR(j)));
grid on;

% Adjust figure size
fig = gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)

print(sprintf('Result/plot/Orthogonal/Noise_variation/DD_Solv_Ext_orth_4_SNR_%s.png',num2str(SNR(j))),'-dpng','-r300')


end


%% === USE CASE III: Skewed noise distribution across wavenumber==========================
% (more noise towards 400cm-1 range and less towards higher wavenumber 1500cm-1 range)


jj=[0.1,0.5,1.0];S_Out2={};SNR=[];

for j =1:length(jj)%[0.01,0.1,0.5,1.5,2.0]
% Add Gaussian noise to the combined signal
noise_std = jj(j)*mean(std(combined_data));%0.1; % Standard deviation for Gaussian noise

factr_noise=[2,5,10];

for kk=1:length(factr_noise)

    noise_std_var = noise_std*ones(length(lam_ftir),1);
    indices2 = find(lam_ftir <= 1500);
    noise_std_var(indices2) = noise_std*linspace(factr_noise(kk),1,length(indices2));

    combined_y = combined_data + noise_std_var .* randn(size(combined_data));

    SNR(j,kk) = mean(mean(combined_y)./std(combined_y));

    Z_ftir=NaN(8,42,1764); 
    Data_ftir = combined_y;

    %Z_ftir(1,1,:)=Data_ftir(:,1)';
    Z_ftir(2,2:9,:)=Data_ftir(:,2:9)';
    Z_ftir(3,10:15,:)=Data_ftir(:,10:15)';
    Z_ftir(4,16,:)=Data_ftir(:,16)';
    Z_ftir(5,17:23,:)=Data_ftir(:,17:23)';
    Z_ftir(6,24:28,:)=Data_ftir(:,24:28)';
    Z_ftir(7,29:35,:)=Data_ftir(:,29:35)';
    Z_ftir(8,36:end,:)=Data_ftir(:,36:end)';
    Z_f=Z_ftir;
    Z_ftir(isnan(Z_ftir))=0;

    P_ftir=(Z_ftir~=0);
    P_ftir=tensor(P_ftir);
    
    Z_ftir=tensor(Z_ftir);
    % figure();plot(noise_std_var)



    SSE_ftir_s2=[];M_FTIR_s2 = struct('tensrs',[]);
    
    % %solvent information
    % %if solvent spectra is known
    Ad = D_ftir_water_sampled;

    S_O_2={};A_O={};

    %_________________________first solvent extraction, later deconvolution____
% for i=1:50
%     opt_options    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',10);
% 
%     R2 = 2;
%     [M_ftir1,~,S_O_2{i},output_ftir] = cp_wopt_s22(Z_ftir, P_ftir, R2,Ad,'lower',0,'opt_options',opt_options);%,'init', M_ftir.u);
% 
%     Z_ftir_s = Z_ftir -tensor(S_O_2{i});
% 
%     opts    = struct( 'factr', 1e0, ...
%     'pgtol', 1e-9, ...
%     'm', 10, ...
%     'maxIts',10000, ...
%     'maxTotalIts',50000, ...
%     'printEvery',500);
% 
%     [M_ftir,~,output_ftir2] = cp_wopt(Z_ftir_s, P_ftir, R_f,'lower',0,'opt_options',opts);%,'init', Minit_ftir);  
% 
%     % [M_ftir,~,S_O_2{i},output_ftir] = cp_wopt_s2(Z_ftir, P_ftir, R_f,j*D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
% %     [M_ftir,~,S_O_2{i},A_O{i},output_ftir] = cp_wopt_s3(Z_ftir, P_ftir, R_f,D_ftir_water_sampled,'lower',0,'init', Minit_ftir);
%     SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(tensor(M_ftir)-P_ftir.*Z_ftir_s)).^2)))];
%     M_FTIR_s2(i).tensrs=M_ftir;
%     disp(['iteration, ' num2str(i)])
% end
%_____________________________________________________________________________

    %________________________________________________________________________________________________
    %___________________________Orthogonal deconvolution with solvent as additional factor________
    % ____________________________________________________________________________________________
    disp(['Iteration ',num2str([j,kk])])

    for i=1:10

    opts    = struct( 'factr', 1e-3, ...
        'pgtol', 1e-5, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',0);
   
    [M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0,'opt_options',opts,'verbosity',0);%,'init', Minit_ftir);  

    
    % 
    % [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
    opts    = struct( 'factr', 1e-3, ...
        'pgtol', 0e-5, ...
        'm', 5, ...
        'maxIts',50000, ...
        'maxTotalIts',50000, ...
        'printEvery',0);

    M_ftir_init = M_ftir.u;
    for ij=1:size(M_ftir_init,1)
        if ij==3
            M_ftir_init{ij} = horzcat(M_ftir_init{ij},abs(j*Ad)); 
        else
            M_ftir_init{ij} = horzcat(M_ftir_init{ij},abs(matrandnorm(size(M_ftir_init{ij},1),1))); 
        
        end
    end

    % Get solvent Free Spctra
    [M_ftir,~,mu_l_o,output_ftir_s] = cp_wopt_s26(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts,'init', M_ftir_init,'verbosity',0);
        
    % Fix spectra and recalibrate the concentration
    A_s = M_ftir.u{3};
    R_c = R_f+1;
    [M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}},'verbosity',0);

    % figure();plot(output_ftir_c.OptOut.err); grid on;
    % sum(sum(sum((double(P_ftir.*tensor(M_ftir)-Z_ftir)).^2)))

    SSE_ftir_s2=[SSE_ftir_s2;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];
    M_FTIR_s2(i).tensrs=M_ftir_c;
%     disp(['iteration, ' num2str(i)])
%     disp(output_ftir_c)
    
    end
    %________________________________________________________________

    % % % csvwrite('Result/SSE_ftir_s.csv',SSE_ftir_s);
    disp('Data deconvolution for intial guess')
    disp(output_ftir)
    disp('Data deconvolution after solvent extraction')
    disp(output_ftir_s)
    disp('Data deconvolution after calibration')
    disp(output_ftir_c)
    
    [sse_minf,sse_locf]=min(SSE_ftir_s2);

    M_ftir_s2=M_FTIR_s2(sse_locf).tensrs;
    % S_Out2{j} = S_O_2{sse_locf};

    % error analysis
    M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s2.u{1}(:,1:R_f);M_ftir_s2.u{2}(:,1:R_f);M_ftir_s2.u{3}(:,1:R_f)});
    rmse_s2 = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
    T_A = reshape(M_ftir_s2.u{3}(:,1:R_f),[size(M_ftir_s2.u{3},1),R_f]);
    % rmse_s2 = sqrt(min(SSE_ftir_s2));
    % T_A = reshape(M_ftir_s2.u{3}(:,1:R_f),[1764,4]);
    profile1 = D_ftir_water_sampled; 
    profile2 = T_A;%mean(mean(double(T_A)));
    Mcos_l_s2 = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
%     disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_s2)])
%     disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_s2)])


    %convert to a ttensor
    R1 = length(M_ftir_s2.lambda);  %<-- Number of factors in X.
    core1 = tendiag(M_ftir_s2.lambda, repmat(R1,1,ndims(M_ftir_s2))); %<-- Create a diagonal core.
    Y_ftir_s2 = ttensor(core1, M_ftir_s2.U) ;%<-- Assemble the ttensor.


    figure()
    ax = gca;
    pos = ax.Position;
    x_pos = pos(1) + 0.1;
    y_pos = pos(2) + 0.6;

    subplot(3,1,1)
    plot(lam_ftir,5*Y_ftir_s2.u{3}(:,1:R_f),lam_ftir, D_ftir_water_sampled,'r*');
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra(SNR= %.2f,noise factor=%.2f)",SNR(j,kk),factr_noise(kk)));
    annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_s2,Mcos_l_s2]), ...
        'FitBoxToText', 'on');
    grid on;

    subplot(3,1,2)
    % Y2 = tensor(S_Out2);%P_ftir.*full(ktensor(S_Out));
    Y2= P_ftir.*ktensor({M_ftir_s2.u{1}(:,R_f+1);M_ftir_s2.u{2}(:,R_f+1);M_ftir_s2.u{3}(:,R_f+1)});
    % S_Out = {Y_ftir_s2.u{1}(:,R_f+1);Y_ftir_s2.u{2}(:,R_f+1);Y_ftir_s2.u{3}(:,R_f+1)};
    % Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;
    mean_tensor = zeros(size(Y2,3),1);
    for i=1:size(Y2,3)
        Y_ext = double(Y2(:,:,i));
           Y_ext(Y_ext == 0) = NaN;
        mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
    end
    plot(1:length(D_ftir_water_sampled),[D_ftir_water_sampled,mean_tensor])
    legend('Ground truth','Mean_Extracted','Interpreter', 'none','Location', 'best')
    title(sprintf("Solvent Substraction,SNR= %.2f",SNR(j,kk)));
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,noise_std_var .* randn(size(combined_data)))
    title('Varying noise level across wavenumber')
    grid on;
    
    % Adjust figure size
    fig = gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)

    print(sprintf('Result/plot/Orthogonal/Noise_variation/factor_var/DD_Solv_Ext_orth_SNR_%.2f_factor_var_%s.png',SNR(j,kk),num2str(factr_noise(kk))),'-dpng','-r300')


end


end


%% ====USE CASE IV: substract solvent effect from a specified region====================================

SSE_ftir=[]; M_FTIR_s = struct('tensrs',[]);

% solvent information--> Partial profile encoding
Ad = zeros(size(lam_ftir));
index_pt = find(lam_ftir>2800 & lam_ftir<3800);    
Ad(index_pt) = D_ftir_water_sampled(index_pt);%==================Signature spectra (Original)
% figure();plot(lam_ftir,Ad)

%_________________________________________________________________

%_________________________________________________________________
% Orthogonal CP decomposition
%_________________________________________________________________

for i=1:10

opts    = struct( 'factr', 1e-3, ...
        'pgtol', 1e-5, ...
        'm', 5, ...
        'maxIts',10000, ...
        'maxTotalIts',50000, ...
        'printEvery',500);
   
[M_ftir,~,output_ftir_s] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0,'opt_options',opts,'verbosity',0);%,'init', Minit_ftir);  
% 
% [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
opts    = struct( 'factr', 1e-3, ...
    'pgtol', 0e-5, ...
    'm', 5, ...
    'maxIts',50000, ...
    'maxTotalIts',50000, ...
    'printEvery',100);

M_ftir_init = M_ftir.u;
for j=1:size(M_ftir_init,1)
    if j==3
        M_ftir_init{j} = horzcat(M_ftir_init{j},abs(Ad)); 
    else
        M_ftir_init{j} = horzcat(M_ftir_init{j},abs(matrandnorm(size(M_ftir_init{j},1),1))); 
        
    end
end

% solvent substraction
[M_ftir,~,mu_l_o,output_ftir_c] = cp_wopt_s26(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts,'init', M_ftir_init,'verbosity',0);

% Fix spectra and recalibrate the concentration
A_s = M_ftir.u{3};
R_c = R_f+1;
[M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}},'verbosity',0);

SSE_ftir=[SSE_ftir;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];
M_FTIR_s(i).tensrs=M_ftir_c;
disp(['iteration, ' num2str(i)])
disp(output_ftir_s)
end
% 

% figure();
% plot([output_ftir.OptOut.err]);
% grid on;
% legend('func','grad')


%_________________________________________________________________

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir_s=M_FTIR_s(sse_locf).tensrs;
% S_Out = S_O{sse_locf};


%convert to a ttensor
R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.


% error analysis
M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s.u{1}(:,1:R_f);M_ftir_s.u{2}(:,1:R_f);M_ftir_s.u{3}(:,1:R_f)});
rmse_s2 = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
T_A = reshape(M_ftir_s.u{3}(:,1:R_f),[size(M_ftir_s.u{3},1),R_f]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_s2)])
disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])


% Plotting
%FTIR%

for i=1:R1
    figure()
    subplot(3,1,1)
    plot(Y_ftir_s.u{1}(:,i));
    legend('conc_T','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,2)
    plot(Y_ftir_s.u{2}(:,i)); hold on
    legend('conc_time','Interpreter', 'none','Location', 'best');
    grid on;

    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_s.u{3}(:,i))
    legend('spectra','Interpreter', 'none','Location', 'best');
    sgtitle(sprintf("FTIR data deconvolution for factor:# %d",i));
    grid on;

end


% load 'Result/workspace_DD_solventFreedata_2024-02-26_12-30'
load 'Result/workspace_DD_solventFreedata_R_3_2024-03-06_15-07'

figure()
org_loc = [1,3,2];
for i=1:R_f
    subplot(R_f,1,i)
    plot(lam_ftir,[Y_ftir_s.u{3}(:,i),Y_ftir.u{3}(:,org_loc(i))]);
    legend(sprintf('PC_%s',num2str(i)),'Solv_free_PC','interpreter','none','Location', 'best');
    grid on;
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
  
end



%Saving
fig = gcf;%figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Orthogonal/Partial/DD_solv_ext_PCs_R_3.png','-dpng','-r300')



figure()
subplot(3,1,1)
plot(lam_ftir,M_ftir_init{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
title(sprintf("(Initialisation)"));
grid on;

subplot(3,1,2)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(3,1,3)
plot(lam_ftir,Y_ftir_s.u{3}(:,1:R_f));
legend('PC_1','PC_2','PC_3','Interpreter', 'none');
sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
grid on;


%Saving
fig = gcf;%figure(2);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Orthogonal/Partial/DD_solv_ext_All_R_3.png','-dpng','-r300')





figure()
subplot(2,1,1)
plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, D_ftir_water_sampled,'r*');
legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
% sgtitle(sprintf("(with water)"));
grid on;

subplot(2,1,2)
% Y2 = tensor(S_Out);%P_ftir.*full(ktensor(S_Out));
S_Out = {Y_ftir_s.u{1}(:,R_f+1);Y_ftir_s.u{2}(:,R_f+1);Y_ftir_s.u{3}(:,R_f+1)};
Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;

mean_tensor = zeros(size(Y2,3),1);
for i=1:size(Y2,3)
    Y_ext = double(Y2(:,:,i));
    Y_ext(Y_ext == 0) = NaN;
    mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
end
plot(1:length(Ad),[Ad,mean_tensor]);% when solvent profile is known
legend('Ground truth','Extracted')
sgtitle(sprintf("Solvent Substraction"));
grid on;
ax = gca;
pos = ax.Position;
x_pos = pos(1) + 0.1;
y_pos = pos(2) + 0.6;
annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similaity:%.2f', [rmse_s2,Mcos_l_ws]), ...
    'FitBoxToText', 'on');



%Saving
fig = gcf;%figure(1);%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/plot/Orthogonal/Partial/DD_solv_ext_first_PCs.png','-dpng','-r300')

% current_timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
% file_name = ['Result/workspace_DD_solEx_', current_timestamp, '.mat'];
% save(file_name, 'M_ftir_s','Y_ftir_s','S_Out','lam_ftir','D_ftir_water_sampled');




%%
% Overall Analysis

% Error Analysis
% 1. RMSE error
rmse_l = sqrt(mean(mean((combined_data-combined_y).^2,2))); % original RMSE across the wavelength dimension for ftir with water spectra data
disp(['RMSE of noise ',num2str(rmse_l)])

% Deconvolution error : without solvent extraction
rmse_ws = sqrt(min(SSE_ftir_ws));
disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])

% Deconvolution error : with solvent extraction
% spectra is known
rmse_sf = sqrt(min(SSE_ftir_s));
disp(['RMSE when deconvoluting data with known solvent spectra extraction ',num2str(rmse_sf)])

% 
% Z = Z_ftir.*P_ftir;
% B = P_ftir.*full(ktensor(M_ftir_ws.u));
% 
% f = 0.5 * norm(Z)^2 - innerprod(Z,B) + 0.5 * norm(B)^2;

% 2. Cosine error (disimilarity)
% original
T_A = tensor(M_ftir);
T_B = Z_ftir;
profile1 = mean(combined_data,2); 
profile2 = mean(combined_y,2);
Mcos_l = 1 - (dot(profile1,profile2) / sqrt(dot(profile1,profile1) * dot(profile2,profile2)));
disp(['Mean cosine error of noise ',num2str(Mcos_l)])

% Data deconvolution without solvent extraction
T_A = reshape(M_ftir_ws.u{3},[1764,4]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_ws = mean(1-((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));

disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])

% Data deconvolution with solvent extraction
% spectra is known
T_A = reshape(M_ftir_s.u{3},[1764,4]);
profile1 = D_ftir_water_sampled; 
profile2 = T_A;%mean(mean(double(T_A)));
Mcos_l_sf = mean(1-((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
disp(['Mean cosine error when deconvoluting data with known solvent spectra extraction ',num2str(Mcos_l_sf)])


% sgtitle("FTIR data deconvolution");
for i=1:R1
    figure()
    plot(lam_ftir,Y_ftir.u{3}(:,i),lam_ftir,Y_ftir_ws.u{3}(:,i),lam_ftir,Y_ftir_s.u{3}(:,i))
    legend('solvent free data','w/o solvent extraction','solvent extraction','Interpreter', 'none','Location', 'best');
    title(sprintf("FTIR data deconvolution for factor:# %d",i)); 
    grid on;
end






