function Fun_SolvSub_V1(Z_ftir, P_ftir,Data_name,lam_ftir,D_ftir_solvent,dir_save, R_f,method,opts)

Ad = 1.0*D_ftir_solvent;%==================Signature spectra (Original)
R2 = 1;
lb=0;
global N_I 


%% ===Data deconvolution Without Solvent extraction=================================
% Fit a PARAFAC model for Z_ftir

WO_solv_substraction(Z_ftir, P_ftir, R_f,Data_name,lam_ftir,D_ftir_solvent,dir_save,lb,opts,method)

%% ===Data deconvolution With Solvent extraction: hard modelling=================================
% solvent prior infromation : full

%solvent information
%if solvent spectra is known
% 
% opt_options    = struct( 'factr', 1e0, ...
% 'pgtol', 1e-9, ...
% 'm', 10, ...
% 'maxIts',10000, ...
% 'maxTotalIts',50000, ...
% 'printEvery',10);


if strcmp(method,'direct') || strcmp(method,'DirectSubstraction') 
    Solv_Sub_dir(Z_ftir, P_ftir, R_f,R2,Data_name,lam_ftir,Ad,dir_save,lb,opts);

    
elseif strcmp(method,'orthogonal')
    Solv_Sub_ortho(Z_ftir, P_ftir, R_f,Data_name,lam_ftir,Ad,dir_save,lb,opts);
    
else
    Disp('wrong method type')
    return;
end
% figure();
% plot([output_ftir.OptOut.err]);
% grid on;
% legend('func','grad')
end

%% Solvent Substraction approaches
function WO_solv_substraction(Z_ftir, P_ftir, R_f,Data_name,lam_ftir,D_ftir_solvent,dir_save,lb,opts,method)
    global N_I
    SSE_ftir_ws=[];M_FTIR_ws = struct('tensrs',[]);
    
    dir_save = strcat(dir_save,'WO_Solv_Subs','/');
    
    if exist(dir_save, 'dir') == 0
        mkdir(dir_save);
        disp(['Directory created: ', dir_save]);
    end

    % opts    = struct( 'factr', 1e-3, ...
    %         'pgtol', 1e-5, ...
    %         'm', 5, ...
    %         'maxIts',10000, ...
    %         'maxTotalIts',50000, ...
    %         'printEvery',500);
    
    for i=1:N_I
        [M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',lb,'opt_options',opts);%,'init', Minit_ftir);
        SSE_ftir_ws=[SSE_ftir_ws;sum(sum(sum((double(P_ftir.*tensor(M_ftir)-Z_ftir)).^2)))];
        M_FTIR_ws(i).tensrs=M_ftir;
        disp(['iteration, ' num2str(i)])
    end
    
    [~,sse_locf]=min(SSE_ftir_ws);
    
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
    T_A = reshape(M_ftir_ws.u{3},[length(lam_ftir),R_f]);
    profile1 = D_ftir_solvent; 
    profile2 = T_A;%mean(mean(double(T_A)));
    Mcos_l_ws = mean(((profile1'*profile2) ./sqrt((profile1'*profile1)*diag(profile2'*profile2))'));
    disp(['RMSE when deconvoluting data without solvent extraction  ',num2str(rmse_ws)])
    disp(['Mean cosine error when deconvoluting data without solvent extraction ',num2str(Mcos_l_ws)])
    
    % Plotting
    %FTIR%
    % 
    % for i=1:R1
    %     figure()
    %     subplot(3,1,1)
    %     plot(Y_ftir_ws.u{1}(:,i));
    %     legend('conc_T','Interpreter', 'none','Location', 'best');
    %     grid on;
    % 
    %     subplot(3,1,2)
    %     plot(Y_ftir_ws.u{2}(:,i));
    %     legend('conc_time','Interpreter', 'none','Location', 'best');
    %     grid on;
    % 
    %     subplot(3,1,3)
    %     plot(lam_ftir,Y_ftir_ws.u{3}(:,i));
    %     legend('spectra','Location', 'best');
    %     grid on;
    %     sgtitle(sprintf("FTIR data deconvolution without solvent extraction for factor:# %d",i));
    %     
    % end
    % 
    
    figure()
    subplot(2,1,1)
    plot(lam_ftir,Y_ftir_ws.u{3});
    legend('PC_1','PC_2','PC_3','PC_4','Interpreter', 'none');
    grid on;

    subplot(2,1,2)
    plot(lam_ftir,Y_ftir_ws.u{3},lam_ftir, D_ftir_solvent,'r*');
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    sgtitle(sprintf("FTIR data deconvolution without solvent extraction: wavenumber vs Spectra"));
    grid on;
    ax = gca;
    pos = ax.Position;
    x_pos = pos(1) + 0.1;
    y_pos = pos(2) + 0.6;
    annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similarity:%.2f', [rmse_ws,Mcos_l_ws]), ...
        'FitBoxToText', 'on');
    
    
    %Saving
    fig = gcf;%gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
    print(strcat(dir_save,strrep(strrep(Data_name,'.mat','.png'),method,'_WOSolvExt_PCs')),'-dpng','-r300')
    
    
    %save data
    file_name = strcat(dir_save,'WO_SolvExt_',Data_name);
    save(file_name, 'M_ftir_ws','Y_ftir_ws','lam_ftir','D_ftir_solvent');
    
    file_name = strcat(dir_save,strrep(Data_name,'.mat',''),'_WOSolvExt_PCs.csv');
    csvwrite(file_name,horzcat(lam_ftir,M_ftir_ws.U{3}));

end

%-----------------------------------------------------------------------------------------
%% ==============================Approach 1: Direct Substraction===============================================
%-----------------------------------------------------------------------------------------
function Solv_Sub_dir(Z_ftir, P_ftir, R_f,R2,Data_name,lam_ftir,Ad,dir_save,lb,opt_options)
    %%_________________________________________________________________
    % First only solvent substraction from the tensor data, later data
    % deconvolution on solvent free tensor data
    %_________________________________________________________________
    % *****1.first solvent substraction and later latent data deconvolution
    global N_I
    SSE_ftir=[];M_FTIR_s = struct('tensrs',[]);
    S_O={};

    dir_save = strcat(dir_save,'DirectSubstraction','/');
    
    if exist(dir_save, 'dir') == 0
        mkdir(dir_save);
        disp(['Directory created: ', dir_save]);
    end


    for i=1:N_I
        [M_ftir1,~,S_O{i},output_ftir1] = cp_wopt_s22(Z_ftir, P_ftir, R2,Ad,'lower',lb,'opt_options',opt_options);%,'init', M_ftir.u);
    
        Z_ftir_s = Z_ftir -tensor(S_O{i});
        % %
        % opts    = struct( 'factr', 1e0, ...
        %     'pgtol', 1e-9, ...
        %     'm', 10, ...
        %     'maxIts',10000, ...
        %     'maxTotalIts',50000, ...
        %     'printEvery',500);
    
        [M_ftir,~,output_ftir] = cp_wopt(Z_ftir_s, P_ftir, R_f,'lower',lb,'opt_options',opt_options);%,'init', Minit_ftir);
    
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
    [~,sse_locf]=min(SSE_ftir);

    M_ftir_s=M_FTIR_s(sse_locf).tensrs;
    S_Out = S_O{sse_locf};
    
    %convert to a ttensor
    R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
    core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
    Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.
    % error analysis
    %-------------Otherwise-------------------------------%
    rmse_s2 = sqrt(min(SSE_ftir));
    T_A = reshape(M_ftir_s.u{3}(:,1:R_f),[size(M_ftir_s.u{3},1),R_f]);
    %-----------------------------------------------------%

    profile1 = Ad; 
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
    
    
    figure()
    subplot(2,1,1)
    plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, Ad,'r*');
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    % sgtitle(sprintf("(with water)"));
    grid on;
    
    subplot(2,1,2)
    Y2 = tensor(S_Out);%P_ftir.*full(ktensor(S_Out));
    
    mean_tensor = zeros(size(Y2,3),1);
    for i=1:size(Y2,3)
        Y_ext = double(Y2(:,:,i));
        Y_ext(Y_ext == 0) = NaN;
        mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
    end
    plot(1:length(Ad),[Ad,mean_tensor]);% when solvent profile is known
    legend('Water Spectra','Extracted')
    sgtitle(sprintf("Solvent Substraction"));
    grid on;
    ax = gca;
    pos = ax.Position;
    x_pos = pos(1) + 0.1;
    y_pos = pos(2) + 0.6;
    annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similaity:%.2f', [rmse_s2,Mcos_l_ws]), ...
        'FitBoxToText', 'on');
    
    
    
    %Saving
    fig = gcf;%gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
    print(strcat(dir_save,Data_name,'_solv_ext_first_PCs.png'),'-dpng','-r300')
    


    figure()
%     subplot(3,1,1)
%     plot(lam_ftir,M_ftir_init{3}(:,1:R_f));
%     legend('PC_1','PC_2','PC_3','Interpreter', 'none');
%     title(sprintf("(Initialisation)"));
%     grid on;
%     
    subplot(2,1,1)
    plot(lam_ftir,Y_ftir_s.u{3},lam_ftir,mean_tensor);
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    % sgtitle(sprintf("(with water)"));
    grid on;
    
    subplot(2,1,2)
    plot(lam_ftir,Y_ftir_s.u{3}(:,1:R_f));
    legend('PC_1','PC_2','PC_3','Interpreter', 'none');
    sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
    grid on;
    
    
    %Saving
    fig = gcf;%gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
    print(strcat(dir_save,Data_name,'_solv_ext_All.png'),'-dpng','-r300')
    
    %save data
    file_name = strcat(dir_save,Data_name);
    save(file_name, 'M_ftir_s','Y_ftir_s','S_Out','lam_ftir','Ad');
    
    file_name = strcat(dir_save,strrep(Data_name,'.mat',''),'.csv');
    csvwrite(file_name,horzcat(lam_ftir,M_ftir_s.U{3},S_Out{3}));


end


% 
% opts    = struct( 'factr', 1e-3, ...
%         'pgtol', 1e-5, ...
%         'm', 5, ...
%         'maxIts',10000, ...
%         'maxTotalIts',50000, ...
%         'printEvery',500);

%-----------------------------------------------------------------------------------------
%% ==============================Approach 2 : Orthogonal Constraints ===============================================
%-----------------------------------------------------------------------------------------

function Solv_Sub_ortho(Z_ftir, P_ftir, R_f,Data_name,lam_ftir,Ad,dir_save,lb,opts)
        
    %_________________________________________________________________
    % Orthogonal CP decomposition
    %_________________________________________________________________
    global N_I
    SSE_ftir=[];M_FTIR_s = struct('tensrs',[]);


    dir_save = strcat(dir_save,'Orthogonal','/');
    
    if exist(dir_save, 'dir') == 0
        mkdir(dir_save);
        disp(['Directory created: ', dir_save]);
    end

    
    for i=1:N_I
       
    [M_ftir,~,output_ftir] = cp_wopt(Z_ftir, P_ftir, R_f,'lower',lb,'opt_options',opts);%,'init', Minit_ftir);  
    % 
    % [M_ftir,~,output_ftir] = cp_wopt_s24(Z_ftir, P_ftir, R_f,Ad,'lower',0,'opt_options',opts);%,'init', M_ftir.u);
%     opts    = struct( 'factr', 1e-3, ...
%         'pgtol', 0e-5, ...
%         'm', 5, ...
%         'maxIts',50000, ...
%         'maxTotalIts',50000, ...
%         'printEvery',100);
    
    M_ftir_init = M_ftir.u;
    for j=1:size(M_ftir_init,1)
        if j==3
            M_ftir_init{j} = horzcat(M_ftir_init{j},abs(Ad)); 
        else
            M_ftir_init{j} = horzcat(M_ftir_init{j},abs(matrandnorm(size(M_ftir_init{j},1),1))); 
            
        end
    end
    
    % Get solvent free spectra
    [M_ftir,~,mu_l_o,output_ftir] = cp_wopt_s26(Z_ftir, P_ftir, R_f,Ad,'lower',lb,'opt_options',opts,'init', M_ftir_init);
    
    % Fix spectra and recalibrate the concentration
    A_s = M_ftir.u{3};
    R_c = R_f+1;
    [M_ftir_c,~,output_ftir_c] = cp_wopt_s27(Z_ftir, P_ftir, R_c,A_s,'lower',0,'opt_options',opts,'init', {M_ftir.u{1},M_ftir.u{2}},'verbosity',0);
    
    SSE_ftir=[SSE_ftir;sum(sum(sum((double(P_ftir.*tensor(M_ftir_c)-Z_ftir)).^2)))];
    
    M_FTIR_s(i).tensrs=M_ftir_c;
    disp(['iteration, ' num2str(i)])
    disp(output_ftir)
    end


    [~,sse_locf]=min(SSE_ftir);

    M_ftir_s=M_FTIR_s(sse_locf).tensrs;
    %convert to a ttensor
    R1 = length(M_ftir_s.lambda);  %<-- Number of factors in X.
    core1 = tendiag(M_ftir_s.lambda, repmat(R1,1,ndims(M_ftir_s))); %<-- Create a diagonal core.
    Y_ftir_s = ttensor(core1, M_ftir_s.U) ;%<-- Assemble the ttensor.

    % error analysis
    %-------------For orthogonal verison------------------%
    M_ftir_ws_PC = P_ftir.*ktensor({M_ftir_s.u{1}(:,1:R_f);M_ftir_s.u{2}(:,1:R_f);M_ftir_s.u{3}(:,1:R_f)});
    rmse_s2 = sqrt(sum(sum(sum((double(M_ftir_ws_PC-Z_ftir)).^2))));%sqrt(min(SSE_ftir_s2));
    T_A = reshape(M_ftir_s.u{3}(:,1:R_f),[size(M_ftir_s.u{3},1),R_f]);
    %-----------------------------------------------------%

    figure()
    subplot(3,1,1)
    plot(lam_ftir,M_ftir_init{3}(:,1:R_f));
    legend('PC_1','PC_2','PC_3','Interpreter', 'none');
    title(sprintf("(Initialisation)"));
    grid on;
    
    subplot(3,1,2)
    plot(lam_ftir,Y_ftir_s.u{3});
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    % sgtitle(sprintf("(with water)"));
    grid on;
    
    subplot(3,1,3)
    plot(lam_ftir,Y_ftir_s.u{3}(:,1:R_f));
    legend('PC_1','PC_2','PC_3','Interpreter', 'none');
    sgtitle(sprintf("FTIR data deconvolution: wavenumber vs Spectra"));
    grid on;
    
    
    %Saving
    fig = gcf;%gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
    print(strcat(dir_save,Data_name,'_solv_ext_All.png'),'-dpng','-r300')


    profile1 = Ad; 
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
    
    
    figure()
    subplot(2,1,1)
    plot(lam_ftir,Y_ftir_s.u{3},lam_ftir, Ad,'r*');
    legend('PC_1','PC_2','PC_3','PC_4','water','Interpreter', 'none');
    % sgtitle(sprintf("(with water)"));
    grid on;
    
    subplot(2,1,2)
    S_Out = {Y_ftir_s.u{1}(:,R_f+1);Y_ftir_s.u{2}(:,R_f+1);Y_ftir_s.u{3}(:,R_f+1)};
    Y2 = P_ftir.*full(ktensor(S_Out));% for orthogonal case;
    
    mean_tensor = zeros(size(Y2,3),1);
    for i=1:size(Y2,3)
        Y_ext = double(Y2(:,:,i));
        Y_ext(Y_ext == 0) = NaN;
        mean_tensor(i) = mean(mean(Y_ext,'omitnan'),'omitnan');
    end
    plot(1:length(Ad),[Ad,mean_tensor]);% when solvent profile is known
    legend('Water Spectra','Extracted')
    sgtitle(sprintf("Solvent Substraction"));
    grid on;
    ax = gca;
    pos = ax.Position;
    x_pos = pos(1) + 0.1;
    y_pos = pos(2) + 0.6;
    annotation('textbox', [x_pos, y_pos, 0.1, 0.1], 'String', sprintf('RMSE: %.2f\n Cosine Similaity:%.2f', [rmse_s2,Mcos_l_ws]), ...
        'FitBoxToText', 'on');
    
    
    
    %Saving
    fig = gcf;%gcf; % Get current figure handle
    fig.PaperUnits = 'inches'; % Set paper units to inches
    fig.PaperPosition = [0 0 13 7]; % Set paper size (8x6 inches, for example)
    print(strcat(dir_save,Data_name,'_solv_ext_first_PCs.png'),'-dpng','-r300')
    
    
    %save data
    file_name = strcat(dir_save,Data_name);
    save(file_name, 'M_ftir_s','Y_ftir_s','S_Out','lam_ftir','Ad');
    
    file_name = strcat(dir_save,strrep(Data_name,'.mat',''),'.csv');
    csvwrite(file_name,horzcat(lam_ftir,M_ftir_s.U{3}));

end