% This script perform following tasks
% > First, based on the WIS and WS, it append all raw data into the these 2 categories
% > Now the data is splitted into 0,10,and 20MeOH categories MANUALLY
% > latet this script preprocess these data files for baseline shift and
% smoothens it further--> rangewise multiple hyperparameter for baseline
% shift applied
%%
clc;clear all; close all;


dir_save='Result/';

if exist(dir_save, 'dir') == 0
    mkdir(dir_save);
    disp(['Directory created: ', dir_save]);
end


if exist('Result/preprocessed/WS', 'dir') == 0
    mkdir('Result/preprocessed/WS');
    disp(['Directory created: ', 'Result/preprocessed/WS']);
end


if exist('Result/preprocessed/WIS', 'dir') == 0
    mkdir('Result/preprocessed/WIS');
    disp(['Directory created: ', 'Result/preprocessed/WIS']);
end

if exist('Result/plot/', 'dir') == 0
    mkdir('Result/plot/');
    disp(['Directory created: ', 'Result/plot/']);
end
%% Append the data
% Specify the directory containing the CSV files
directory = 'K:\SolventSubtraction\BioOil\Data\FTIR data';  % Replace 'path_to_directory' with your directory path

%------------------------for WIS---------------------------------------
% Get a list of all CSV files with a specific prefix

filePattern = fullfile(directory, 'WIS *.csv');  % 
fileList = dir(filePattern);
fileList.name

% Initialize variables to store x and y data
xData = [];
yData = cell(numel(fileList), 1);  % Cell array to store y data from each file

% Loop through each CSV file and extract data

yData_com=[];

for i = 1:numel(fileList)
    filename = fullfile(directory, fileList(i).name);
    data = csvread(filename);
    
    disp(['File Name ',fileList(i).name])
    length(data)
    % Extract x data from the current file
    xData_current = data(:, 1);
    
    % Check if x data is the same for all files
    if isempty(xData) || isequal(xData, xData_current)
        xData = xData_current;  % Update x data if it's the same
        yData{i} = data(:, 2);  % Store y data from the current file

        yData_com = [yData_com,data(:, 2)];
    else
        idx = find(ismember(xData_current, xData));
        
        if isempty(idx)
            yData{i}=NaN*ones(size(yData_com,1),1);
        else
            yData{i}=NaN*ones(size(yData_com,1),1);
            yData{i}(idx) = data(idx,2);
        end
         yData_com = [yData_com,yData{i}];
    end
    
    % Optionally, you can display the data from each file
%     disp(['Data from ' fileList(i).name ':']);
%     disp(data);

%     yData_com = [yData_com,yData{i}(:,1)];
end

% Combine x data and appended y data

BM_WIS_data = [xData, yData_com];

% Display the combined data
disp('Combined Data:');
disp(BM_WIS_data);

% Save Data

columnNames = {'wavenumber',fileList.name};
for i = 1:numel(columnNames)
    columnNames{i} = strrep(columnNames{i}, '.CSV', '');
    columnNames{i} = strrep(columnNames{i}, 'WIS', '');
    
end

% Create a table with the data and column names
dataTable = array2table(BM_WIS_data, 'VariableNames', columnNames);

% Save the table to a CSV file
writetable(dataTable, 'K:\SolventSubtraction\BioOil\Data\FTIR data\Append\WIS_sampled_data.csv');


%------------------------for WS---------------------------------------
% Get a list of all CSV files with a specific prefix
filePattern = fullfile(directory, 'WS *.csv');  % 'WS * 0_*.csv'
fileList = dir(filePattern);

% Initialize variables to store x and y data
xData = [];
yData = cell(numel(fileList), 1);  % Cell array to store y data from each file

% Loop through each CSV file and extract data

yData_com=[];

for i = 1:numel(fileList)
    filename = fullfile(directory, fileList(i).name);
    data = csvread(filename);
    
    disp(['File Name ',fileList(i).name])
    length(data)
    % Extract x data from the current file
    xData_current = data(:, 1);
    
    % Check if x data is the same for all files
    if isempty(xData) || isequal(xData, xData_current)
        xData = xData_current;  % Update x data if it's the same
        yData{i} = data(:, 2);  % Store y data from the current file

        yData_com = [yData_com,data(:, 2)];
    else
        disp(['change in number of data points for file ',fileList(i).names])

        idx = find(ismember(xData_current, xData));
        
        if isempty(idx)
            yData{i}=NaN*ones(size(yData_com,1),1);
        else
            yData{i}=NaN*ones(size(yData_com,1),1);
            yData{i}(idx) = data(idx,2);
        end
         yData_com = [yData_com,yData{i}];
    end
    
    % Optionally, you can display the data from each file
%     disp(['Data from ' fileList(i).name ':']);
%     disp(data);

%     yData_com = [yData_com,yData{i}(:,1)];
end

% Combine x data and appended y data

BM_WS_data = [xData, yData_com];

% Display the combined data
disp('Combined Data:');
disp(BM_WS_data);

% Save Data

columnNames = {'wavenumber',fileList.name};
for i = 1:numel(columnNames)
    columnNames{i} = strrep(columnNames{i}, '.CSV', '');
    columnNames{i} = strrep(columnNames{i}, 'WS', '');
    
end

% Create a table with the data and column names
dataTable = array2table(BM_WS_data, 'VariableNames', columnNames);

% Save the table to a CSV file
writetable(dataTable, 'K:\SolventSubtraction\BioOil\Data\FTIR data\Append\WS_sampled_data.csv');


%% Directories

% -----------------solvent-----------------------------

filename_ftir_water = 'Data/FTIR data/Pure water.CSV';
filename_ftir_MeOH = 'Data/FTIR data/MeOH.CSV';
filename_ftir_10_MeOH = 'Data/FTIR data/10_ MeOH in water.CSV';
filename_ftir_20_MeOH = 'Data/FTIR data/20_ MeOH in water.CSV';

Data_ftir_water= xlsread(filename_ftir_water);  %data matrix
lam_ftir_water = Data_ftir_water(:,1);
D_ftir_water = Data_ftir_water(:,2);

Data_ftir_MeOH= xlsread(filename_ftir_MeOH);  %data matrix
lam_ftir_MeOH = Data_ftir_MeOH(:,1);
D_ftir_MeOH = Data_ftir_MeOH(:,2);

Data_ftir_10MeOH= xlsread(filename_ftir_10_MeOH);  %data matrix
lam_ftir_10MeOH = Data_ftir_10MeOH(:,1);
D_ftir_10MeOH = Data_ftir_10MeOH(:,2);

Data_ftir_20MeOH= xlsread(filename_ftir_20_MeOH);  %data matrix
lam_ftir_20MeOH = Data_ftir_20MeOH(:,1);
D_ftir_20MeOH = Data_ftir_20MeOH(:,2);

figure()
plot(lam_ftir_water,D_ftir_water,'LineWidth', 3);hold on
plot(lam_ftir_MeOH,D_ftir_MeOH,'LineWidth', 3)
plot(lam_ftir_10MeOH,D_ftir_10MeOH,'r*',MarkerSize=3)
plot(lam_ftir_20MeOH,D_ftir_20MeOH,'bd',MarkerSize=3)
grid on;
legend('Water','MeOH','10_MeOH','20_MeOH','Interpreter', 'none')
title('Solvent FTIR spectra')

% Save figure as PDF
fig = gcf;%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print('Result/SolventFTIR.png', '-dpng', '-r300');

%%

%-----------------------------reaction data----------------
tag= 'WIS'; % WIS, WS
files = dir(fullfile(pwd,strcat('\Data\FTIR data\Append\',tag,'\',tag,'*.xlsx')));

for j=1:numel(files)
filename_ftir = strcat('Data\FTIR data\Append\',tag,'\',files(j).name);

disp(['Iteration ',num2str(j),' ','Filename ',files(j).name])

D_ftir= xlsread(filename_ftir);  %data matrix
time_ftir=unique(D_ftir(1,2:end)); %temperature
time_ftir = time_ftir(~isnan(time_ftir));
T_ftir=unique(D_ftir(2,2:end)); %time
X_ftir=D_ftir(3:end,:);%intensity
lam_ftir=X_ftir(:,1);%wavenumber
X11_ftir=X_ftir(:,2:end); %intensity
D_ftir=X11_ftir;%2-log10(X11_ftir);% change to absorbance which is interpreted as conc

% subplot(2,1,1)
% plot(lam_ftir,D_ftir(:,1))
% subplot(2,1,2)
% plot(lam_ftir_water,0.05*D_ftir_water)


% Data preprocessing
D_ftir_c = D_ftir;

allNaNColumns = all(isnan(D_ftir_c));

I_notAllNaN = find(~allNaNColumns);


%baseline and background correction
split = find(lam_ftir< 2000, 1, 'last' );
Dback_ftir = D_ftir_c(:,I_notAllNaN);
D_ftir_MeOH_proc = D_ftir_MeOH;
D_ftir_water_proc = D_ftir_water;
D_ftir_10MeOH_proc =D_ftir_10MeOH;
D_ftir_20MeOH_proc = D_ftir_20MeOH;


if strcmp(tag, 'WS')
    Dback_ftir(1:split,:)=msbackadj(lam_ftir(1:split),D_ftir_c(1:split,I_notAllNaN));%,'WindowSize',500); % Correct baseline of signal--> spline approximation--> noise removal
    Dback_ftir(split+1:end,:)=msbackadj(lam_ftir(split+1:end),D_ftir_c(split+1:end,I_notAllNaN),'WindowSize',500); % Correct baseline of signal--> spline approximation--> noise removal

    % Solvent %

    D_ftir_water_proc(1:split)=msbackadj(lam_ftir_water(1:split),D_ftir_water(1:split));%,'WindowSize',400);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    D_ftir_water_proc(split+1:end)=msbackadj(lam_ftir_water(split+1:end),D_ftir_water(split+1:end),'WindowSize',500);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
   
    D_ftir_10MeOH_proc(1:split)=msbackadj(lam_ftir_10MeOH(1:split),D_ftir_10MeOH(1:split));%,'WindowSize',400);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    D_ftir_10MeOH_proc(split+1:end)=msbackadj(lam_ftir_10MeOH(split+1:end),D_ftir_10MeOH(split+1:end),'WindowSize',500);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
   
    D_ftir_20MeOH_proc(1:split)=msbackadj(lam_ftir_20MeOH(1:split),D_ftir_20MeOH(1:split));%,'WindowSize',400);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    D_ftir_20MeOH_proc(split+1:end)=msbackadj(lam_ftir_20MeOH(split+1:end),D_ftir_20MeOH(split+1:end),'WindowSize',500);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
   
else
    Dback_ftir(1:split,:)=msbackadj(lam_ftir(1:split),D_ftir_c(1:split,I_notAllNaN));%,'WindowSize',400);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    Dback_ftir(split+1:end,:)=msbackadj(lam_ftir(split+1:end),D_ftir_c(split+1:end,I_notAllNaN),'WindowSize',475);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal

    % Solvent %

    D_ftir_MeOH_proc(1:split)=msbackadj(lam_ftir_MeOH(1:split),D_ftir_MeOH(1:split));%,'WindowSize',400);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    D_ftir_MeOH_proc(split+1:end)=msbackadj(lam_ftir_MeOH(split+1:end),D_ftir_MeOH(split+1:end),'WindowSize',475);%,'StepSize',150); % Correct baseline of signal--> spline approximation--> noise removal
    
end

Data_ftir=Dback_ftir;%mssgolay(lam_ftir,Dback_ftir,'DEGREE',2,'Span',5);% Smooth signal with peaks using least-squares polynomial--> Savitzky and Golay filters-->addtive noise removal

Noise_ftir=Dback_ftir-Data_ftir;


% Plotting
%FTIR%
figure()
loc = round(linspace(1,size(Data_ftir,2),size(Data_ftir,2)));
sgtitle(sprintf('Data Pre-processing for: %s',strrep(files(j).name,'.xlsx','')),'Interpreter', 'None')
% labels = {'original','baseline adjusted','smoothening'};
D_org = D_ftir_c(:,I_notAllNaN);
for i=1:length(loc)
    subplot(2,ceil((length(loc)+1)/2),i)
    plot(lam_ftir,D_org(:,loc(i)),'DisplayName','original'); hold on
    plot(lam_ftir,Dback_ftir(:,loc(i)),'r*','MarkerSize', 2,'DisplayName','baseline adjusted'); hold on
    plot(lam_ftir,Data_ftir(:,loc(i)),'bd','MarkerSize', 2,'DisplayName','smoothening')

    % Add legend
    legend('Location', 'best');
    
    % Add title
    title(sprintf("FTIR data Preprocessing:# %d",loc(i)));
    
    % Add axis labels
    xlabel('wavenumber');
    ylabel('Absorbance');
    grid on;
end

subplot(2,ceil((length(loc)+1)/2),i+1)
if strcmp(tag, 'WIS')

    plot(lam_ftir_MeOH,D_ftir_MeOH,'MarkerSize', 2,'DisplayName','Solv_MeOH_org'); hold on
    plot(lam_ftir_MeOH,D_ftir_MeOH_proc,'MarkerSize', 2,'DisplayName','Solv_MeOH_proc')
else
    plot(lam_ftir_water,D_ftir_water,'MarkerSize', 2,'DisplayName','Solv_Water_org');hold on
    plot(lam_ftir_10MeOH,D_ftir_10MeOH,'MarkerSize', 2,'DisplayName','Solv_10%_MeOH_org');hold on
    plot(lam_ftir_20MeOH,D_ftir_20MeOH,'MarkerSize', 2,'DisplayName','Solv_20%_MeOH_org');hold on

    plot(lam_ftir_water,D_ftir_water_proc,'MarkerSize', 2,'DisplayName','Solv_Water_proc');hold on
    plot(lam_ftir_10MeOH,D_ftir_10MeOH_proc,'MarkerSize', 2,'DisplayName','Solv_10%_MeOH_proc');hold on
    plot(lam_ftir_20MeOH,D_ftir_20MeOH_proc,'MarkerSize', 2,'DisplayName','Solv_20%_MeOH_proc')
end
% Add legend
legend('Location', 'best',Interpreter='none');

% Add title
title("FTIR data Preprocessing: Solvent");

% Add axis labels
xlabel('wavenumber');
ylabel('Absorbance');
grid on;


% Save figure as PDF
fig = gcf;%gcf; % Get current figure handle
fig.PaperUnits = 'inches'; % Set paper units to inches
fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
print(strcat('Result/preprocessed/',tag,'/','fig_',strrep(files(j).name,'.xlsx',''),'.png'), '-dpng', '-r300');

Data_ftir_org=nan(size(D_ftir));
Data_ftir_org(:,I_notAllNaN) =Data_ftir;
Data_ftir_org =horzcat(lam_ftir,Data_ftir_org);

if strcmp(tag, 'WIS')
    Data_ftir_org_MeOH =horzcat(lam_ftir_MeOH,D_ftir_MeOH_proc);

    file_name_solv = strcat('Result/preprocessed/','D_solv_MeOH', '.mat');
    save(file_name_solv, 'Data_ftir_org_MeOH');

else
    Data_ftir_org_water =horzcat(lam_ftir_water,D_ftir_water_proc);

    file_name_solv = strcat('Result/preprocessed/','D_solv_water', '.mat');
    save(file_name_solv, 'Data_ftir_org_water');

    Data_ftir_org_10MeOH =horzcat(lam_ftir_10MeOH,D_ftir_10MeOH_proc);
    
    file_name_solv = strcat('Result/preprocessed/','D_solv_10MeOH', '.mat');
    save(file_name_solv, 'Data_ftir_org_10MeOH');
    
    Data_ftir_org_20MeOH =horzcat(lam_ftir_20MeOH,D_ftir_20MeOH_proc);

    file_name_solv = strcat('Result/preprocessed/','D_solv_20MeOH', '.mat');
    save(file_name_solv, 'Data_ftir_org_20MeOH');
end



file_name = strcat('Result/preprocessed/',tag,'/','Preprocessed_',strrep(files(j).name,'.xlsx',''), '.mat');
save(file_name, 'Data_ftir_org');

%stadard normal variate

end




