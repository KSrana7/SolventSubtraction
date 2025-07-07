% This script perform following tasks
% Process data files for baseline shift and
% smoothens it further--> rangewise multiple hyperparameter for baseline
% shift applied
%%

clc; clear all; close all;

set(groot,'defaultLineLineWidth',2)

addpath('K:\SolventSubtraction\ToolBoxes\tensor_toolbox-master\tensor_toolbox-master');
addpath('K:\SolventSubtraction\ToolBoxes\L-BFGS-B-C-master\L-BFGS-B-C-master\Matlab');
addpath('K:\SolventSubtraction\ToolBoxes\poblano_toolbox_1.1\poblano_toolbox_1.1');

addpath('K:\SolventSubtraction\ToolBoxes\MBtoolbox_v_02')
addpath('K:\SolventSubtraction\ToolBoxes\nway331')


dir_save='Result/';

if exist(dir_save, 'dir') == 0
    mkdir(dir_save);
    disp(['Directory created: ', dir_save]);
end

%%
% ----------------Data----------------------------

filename_ftir = 'Data/Biomass_combined_ftir.xlsx';
D_ftir= readmatrix(filename_ftir);  %data matrix
T_ftir=D_ftir(1,2:end); 
T_ftir = T_ftir(~isnan(T_ftir));%temperature
time_ftir=D_ftir(2,2:end); 
time_ftir=time_ftir(~isnan(time_ftir));%time
X_ftir=D_ftir(3:end,~isnan(D_ftir(1,:)));%intensity
sum(isnan(X_ftir))

lam_ftir=D_ftir(3:end,1);%wavenumber
% X11_ftir=X_ftir(:,2:end); %intensity

max_v = log10(max(X_ftir));
X_ftir = max_v-log10(X_ftir); % change to absorbance which is interpreted as conc
% X_ftir = 2-log10(X_ftir); % change to absorbance which is interpreted as conc
figure();plot(lam_ftir,X_ftir);grid on

% Remove uninformative bands
split_noise=850;
i_noise = find(lam_ftir<split_noise,1,"last")

lam_ftir=lam_ftir(i_noise:end,1);%wavenumber
X_ftir=X_ftir(i_noise:end,:); %intensity
figure();plot(lam_ftir,X_ftir);grid on

disp('Count of imaginary number');sum(imag(X_ftir)~=0)
% X_ftir = X_ftir(128:end,:); % to ignore the imaginary numbers (not possible)
% lam_ftir = lam_ftir(128:end); % to ignore the imaginary numbers (not possible)
% disp('Count of imaginary number after removal');sum(imag(X_ftir)~=0)

%----Solvent----%
filename_ftir_water = 'Data/FTIR_Water_nist.xlsx';

Data_ftir_water= readmatrix(filename_ftir_water);  %data matrix
lam_ftir_water = Data_ftir_water(:,1);
D_ftir_water = Data_ftir_water(:,4);


disp('Count of imaginary number');sum(imag(D_ftir_water)~=0)


sum(isnan(D_ftir_water))


%% Data preprocessing
% for jj=[200,250,300,350,400,450,500]


D_ftir_c = X_ftir;
Dback_ftir = X_ftir;
split = find(lam_ftir< 2500, 1, 'last' );

series1 = msbackadj(lam_ftir,D_ftir_c,'WindowSize',800);
% series1(:,[4,5,3,8]) = msbackadj(lam_ftir,D_ftir_c(:,[4,5,3,8]),'WindowSize',200);
% series1(:,[7,9]) = msbackadj(lam_ftir,D_ftir_c(:,[7,9]),'WindowSize',100);

series2 = msbackadj(lam_ftir,D_ftir_c,'WindowSize',200);

% series3 = series2;
% series3(:,[5]) = msbackadj(lam_ftir,D_ftir_c(:,5),'WindowSize',170);
% series3(:,[8,3]) = msbackadj(lam_ftir,D_ftir_c(:,[8,3]),'WindowSize',130);
% series3(:,[7]) = msbackadj(lam_ftir,D_ftir_c(:,7),'WindowSize',110);
% series3(:,[9]) = msbackadj(lam_ftir,D_ftir_c(:,9),'WindowSize',110);
% split2 = find(lam_ftir< 997, 1, 'last' );

Dback_ftir(1:split,:)=series2(1:split,:); % Correct baseline of signal--> spline approximation--> noise removal

Dback_ftir(split+1:end,:)=series1(split+1:end,:);% msbackadj(lam_ftir(split+1:end,:),D_ftir_c(split+1:end,:),'WindowSize',ws); % Correct baseline of signal--> spline approximation--> noise removal
% Dback_ftir(1:split2,:)=series3(1:split2,:);% msbackadj(lam_ftir(split+1:end,:),D_ftir_c(split+1:end,:),'WindowSize',ws); % Correct baseline of signal--> spline approximation--> noise removal


% figure();plot(lam_ftir,[Dback_ftir(:,[3,5,7,8,9])]); grid on;legend(string([3,5,7,8,9])) %(string([1,2,6,4,5,3,7,8,9]))
% figure();plot(lam_ftir,[Dback_ftir(:,[1,2,4,6])]); grid on;legend(string([1,2,4,6])) %(string([1,2,6,4,5,3,7,8,9]))


Data_ftir=mssgolay(lam_ftir,Dback_ftir,'DEGREE',2,'Span',5);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal
% Data_ftir=Dback_ftir;
Noise_ftir=Dback_ftir-Data_ftir;

figure();plot(lam_ftir,Data_ftir(:,1:9)); grid on;legend%(string([1,2,3,4,5,6,7,8,9]))

% Plotting
%FTIR%
figure()
loc = round(linspace(1,27,6));
% sgtitle(sprintf("FTIR data Preprocessing for span:# %d",jj))
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
% end
% Save figure as PDF
fig = gcf;%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/Data_preprocess.png','-dpng','-r300')


Dback_ftir_s = D_ftir_water;

% Data_ftir_s(1:split,:)=mssgolay(lam_ftir_water2(1:split),Dback_ftir_s(1:split,:), ...
%     'DEGREE',2,'Span',20);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal
% Data_ftir_s(split+1:end,:)=mssgolay(lam_ftir_water2(split+1:end),Dback_ftir_s(split+1:end,:), ...
%     'DEGREE',2,'Span',100);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal

% Dback_ftir_s(1:split,:)=msbackadj(lam_ftir_water(1:split),Dback_ftir_s(1:split,:),'WindowSize',200); % Correct baseline of signal--> spline approximation--> noise removal

Dback_ftir_s(1:end,:)=msbackadj(lam_ftir_water(1:end),Dback_ftir_s(1:end,:),'WindowSize',1000); % Correct baseline of signal--> spline approximation--> noise removal

figure();plot(lam_ftir_water,Dback_ftir_s,lam_ftir_water,D_ftir_water)


% Dback_ftir_s = D_ftir_water;
Data_ftir_s=mssgolay(lam_ftir_water,Dback_ftir_s,'DEGREE',2,'Span',5);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal

sum(isnan(Data_ftir_s))
% plot(lam_ftir_water,Data_ftir_s,lam_ftir_water,D_ftir_water)

% plot(lam_ftir_water,Data_ftir_s)
% 
% Data_ftir_s=mssgolay(lam_ftir_water2,Dback_ftir_s, ...
%     'DEGREE',2,'Span',100);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal




%% Sample signals

% Tolerance for accepting nearest x-values
tolerance = 0.5;

% Initialize combined signal
combined_lam = [];
combined_data = Data_ftir;
D_ftir_water_sampled=zeros(length(lam_ftir),1);
total = [];

% Find the nearest x-values in signal 2 for each x-value in signal 1
for i = 1:length(lam_ftir)
    % Find the index of the nearest x-value
    [~, idx] = min(abs(lam_ftir_water - lam_ftir(i)));
    
    % Check if the nearest x-value is within tolerance
    if (abs(lam_ftir_water(idx) - lam_ftir(i))) <= tolerance
        combined_lam = [combined_lam, lam_ftir(i)];
        D_ftir_water_sampled(i) = 0.1*D_ftir_water(idx);
%         combined_data(i,:) = combined_data(i,:) + D_ftir_water_sampled(i); % Add the corresponding y-values as per their weight fraction
    else
        % Print exception message
        total = [total, i];
        disp(['No matching x-value within tolerance for x = ', num2str(lam_ftir(i))]);
    end
end
%aaa = combined_data - D_ftir %validation
% Add Gaussian noise to the combined signal
noise_std = 0.0*mean(std(combined_data));%0.1; % Standard deviation for Gaussian noise
combined_y = combined_data + noise_std * randn(size(combined_data));

norm(combined_y-combined_data)

%% Data collection
D_ftir; % original ftir data
D_ftir_water_sampled; % water spectra downsampled to match the wavenumbeer
combined_data;% data with water spectra
combined_y; % data with water spectra and random noise 

% Plot the original and combined signals
figure;
loc = round(linspace(1,27,6));
for i=1:length(loc)
    subplot(2,3,i)
    plot(lam_ftir, D_ftir_water_sampled, 'r*',lam_ftir, Data_ftir(:,i), 'b-');
    legend('water spectra', 'mixture spectra', 'Location', 'best');
    xlabel('waveumber cm-1');
    ylabel('absorbance');
    title(sprintf('Combining Signal at (%d min, %d 0C)', [time_ftir(loc(i)),T_ftir(loc(i))]));
    grid on
end
% Save figure as PDF
fig = gcf;%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/Data_plot_with_water.png','-dpng','-r300')


%% Saving Data


file_name = ['Result/workspace_Data_Pre_Process_', 'V0', '.mat'];
save(file_name,'Data_ftir','D_ftir_water_sampled','lam_ftir','T_ftir','time_ftir','X_ftir','D_ftir_water');

disp('Saved Successfully')