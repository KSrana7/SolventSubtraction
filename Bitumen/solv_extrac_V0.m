% script for solvent extraction: mainly content from tensor_df_cc.m.
% it contains:
%   Data generation that contains fraction of water spectra
%   Data preprocessing--> baseline correction and smoothening
%   Rank determination
%   Data deconvolution--> CP_wopt original and modified version for solvent
%   hard modelling


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
        combined_data(i,:) = combined_data(i,:) + 0.05*D_ftir_water(idx); % Add the corresponding y-values as per their weight fraction
    else
        % Print exception message
        total = [total, i];
        disp(['No matching x-value within tolerance for x = ', num2str(lam_ftir(i))]);
    end
end
%aaa = combined_data - D_ftir %validation
% Add Gaussian noise to the combined signal
noise_std = 0.0*mean(var(combined_data));%0.1; % Standard deviation for Gaussian noise
combined_y = combined_data + noise_std * randn(size(combined_data));

% Plot the original and combined signals
figure;
loc = round(linspace(1,42,6));
for i=1:length(loc)
    subplot(2,3,i)
    plot(lam_ftir_water, 0.05*D_ftir_water, 'r--'); hold on
    plot(lam_ftir, D_ftir(:,i), 'b-');hold on
    plot(lam_ftir, combined_y(:,i), 'k-');
    legend('water spectra', 'mixture spectra', 'Combined spectra', 'Location', 'best');
    xlabel('waveumber cm-1');
    ylabel('absorbance');
    title(sprintf('Combining Signal at (%d min, %d 0C)', [time_ftir(loc(i)),T_ftir(loc(i))]));
    grid on
end
% Save figure as PDF
% print('Result/plot/combined_figure.png', '-dpng', '-r300');

%% Data preprocessing
D_ftir_c = combined_y;

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
% 
% %NMR%
% figure(2)
% loc = round(linspace(1,32,6));
% % labels = {'original','baseline adjusted','smoothening'};
% 
% for i=1:length(loc)
%     subplot(2,3,i)
%     plot(conc_hnmr,D_hnmr(:,loc(i)),'DisplayName','original'); hold on
%     plot(conc_hnmr,Dback_hnmr(:,loc(i)),'DisplayName','baseline adjusted');hold on
%     plot(conc_hnmr,Data_hnmr(:,loc(i)),'DisplayName','smoothening')
%     
%     % Add legend
%     legend('Location', 'best');
%     
%     % Add title
%     title(sprintf("NMR data Preprocessing:# %d",loc(i)));
%     
%     % Add axis labels
%     xlabel('Concentration');
%     ylabel('Intensity');
%     grid on;
% end

% Save figure as PDF
% print('Result/plot/preprocessed_figure.png', '-dpng', '-r300');


%stadard normal variate




%% Determine Rank
global Z_ftir; global P_ftir;
Z_ftir=NaN(8,42,1764); 

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

% print('Result/plot/rank_deter.png','-dpng','-r300')

% [A_f,B_f,C_f]=fac2let(Factors_ftir);
% M_f = nmodel(Factors_ftir);

R_f=4;



%% Data deconvolution
% Fit a PARAFAC model for Z_ftir

% Set up optimization parameters
% Get the defaults
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
% Tighten relative change in function value tolearnce. This is often too large.
ncg_opts.RelFuncTol = 1.0e-20;
% Increase the number of iterations.
ncg_opts.MaxIters = 10^4;
% Only display every 10th iteration
ncg_opts.DisplayIters = 10;
% Display the final set of options
ncg_opts

% Initial guess
%[Ainit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,1)*tenmat(Z_ftir,1)'),R_f,[]);
%[Binit_ftir,~]=NNDSVD(double(tenmat(Z_ftir,2)*tenmat(Z_ftir,2)'),R_f,[]);
%Cinit_ftir=double(tenmat(Z_ftir,3)*khatrirao(Binit_ftir,Ainit_ftir))*pinv((Binit_ftir'*Binit_ftir).*(Ainit_ftir'*Ainit_ftir));
%Minit_ftir={Ainit_ftir;Binit_ftir;Cinit_ftir};
%Minit_ftir = create_guess('Data',Z_ftir, 'Num_Factors', R_f, ...
  %    'Factor_Generator', 'nvecs');
%Call optimizer

SSE_ftir=[];M_FTIR = struct('tensrs',[]);

for i=1:50
[M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',0);%,'init', Minit_ftir);
SSE_ftir=[SSE_ftir;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

csvwrite('Result/SSE_ftir.csv',SSE_ftir);

[sse_minf,sse_locf]=min(SSE_ftir)

M_ftir=M_FTIR(sse_locf).tensrs;


csvwrite('Result/M_ftirU1.csv',M_ftir.U{1});
csvwrite('Result/M_ftirU2.csv',M_ftir.U{2});
csvwrite('Result/M_ftirU3.csv',M_ftir.U{3});

%convert to a ttensor
R1 = length(M_ftir.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir.lambda, repmat(R1,1,ndims(M_ftir))); %<-- Create a diagonal core.
Y_ftir = ttensor(core1, M_ftir.U) ;%<-- Assemble the ttensor.

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
    sgtitle(sprintf("FTIR data deconvolution for factor:# %d",i));
    
end


figure()
plot(Y_ftir.u{1});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none','Location', 'best');
title(sprintf("FTIR data deconvolution: Conc vs Temp"));
grid on;

figure()
plot(Y_ftir.u{2});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none','Location', 'best');
title(sprintf("FTIR data deconvolution: Conc vs time"));
grid on;

figure()
plot(lam_ftir,Y_ftir.u{3});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none','Location', 'best');
title(sprintf("FTIR data deconvolution: wavenumber vs Spctra"));
grid on;



print('Result/plot/parafac.png','-dpng','-r300')
%% Data deconvolution with solvent extraction--> hard modelling of the
%% Data deconvolution
% Fit a PARAFAC model for Z_ftir

% Set up optimization parameters
% Get the defaults
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
% Tighten relative change in function value tolearnce. This is often too large.
ncg_opts.RelFuncTol = 1.0e-20;
% Increase the number of iterations.
ncg_opts.MaxIters = 10^4;
% Only display every 10th iteration
ncg_opts.DisplayIters = 10;
% Display the final set of options
ncg_opts

% solvent prior infromation
SSE_ftir_s=[];M_FTIR_s = struct('tensrs',[]);

%solvent information
Ad = zeros(length(lam_ftir),1);
indices = find(lam_ftir >= 1580 & lam_ftir <= 1680);
Ad(indices) = 1;
indices = find(lam_ftir >= 3200 & lam_ftir <= 3500);
Ad(indices) = 1;

%plot(lam_ftir,Ad,'r--',lam_ftir,D_ftir_water_sampled)
Ad = D_ftir_water_sampled;

for i=1:50
[M_ftir,~,output_ftir] = cp_wopt_s(Z_ftir, P_ftir, R_f+1,Ad,'lower',0);%,'init', Minit_ftir);
SSE_ftir_s=[SSE_ftir_s;sum(sum(sum((double(tensor(M_ftir)-Z_ftir)).^2)))];
M_FTIR_s(i).tensrs=M_ftir;
disp(['iteration, ' num2str(i)])
end

% % % csvwrite('Result/SSE_ftir_s.csv',SSE_ftir_s);

[sse_minf,sse_locf]=min(SSE_ftir_s)

M_ftir_s=M_FTIR_s(sse_locf).tensrs;


% csvwrite('Result/M_ftirU1_s.csv',M_ftir_s.U{1});
% csvwrite('Result/M_ftirU2_s.csv',M_ftir_s.U{2});
% csvwrite('Result/M_ftirU3_s.csv',M_ftir_s.U{3});

%convert to a ttensor
R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.

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

figure()
plot(lam_ftir,Y_ftir_s.u{3});
legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none','Location', 'best');
title(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
grid on;


%print('Result/plot/parafac.png','-dpng','-r300')





figure()
subplot(2,1,1)
plot(lam_ftir,D_ftir.u{3}(1),'r-',lam_ftir,Y_ftir_s.u{3}(1),'b*')





