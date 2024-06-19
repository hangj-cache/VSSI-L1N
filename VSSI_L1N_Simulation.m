clc
clear
tic
%brainstorm %norgui
%%
algorithms = {'VSSI-L1N','SISSY','VB-SCCD','wMNE','LORETA'};
ResultsLoading = [0 0 0 0 0];
% ResultsLoading = [1 1 1 1 1 1 1];
WGN = 1; % Using white Gaussian Nosie (WGN = 1) or Human EEG noise (WGN  = 0);
LFNormlization = 0; % Whether normalize the LeadField Matrix
Uniform = 1; % Uniform/NonUniform Sources
VariousExtents = 0;
VariousSNRs = 0;
VariousSNIRs = 1;
VariousPatches = 0;
VariousCorrelation = 0;
VariousChannels = 0;
Test = 0;
if VariousExtents+VariousSNRs+VariousSNIRs+VariousPatches+VariousCorrelation+VariousChannels+Test ~= 1
    error('There will be one and only one scenario.');
end
% tic
%% Export the Channel,Cortex and HeadModel to workspace
if WGN
    channelselect=[1:32,34:42,44:64]; %Select EEG data
%     channelselect=[1:32 34:42 44:59 61:63];
else
    channelselect=[1:32 34:42 44:59 61:63]; % for Real noise simulations
end
% [sStudy, iStudy] = bst_get('StudyWithCondition', 'LiDaoli/Stim_2');
[sStudy, iStudy] = bst_get('StudyWithCondition','Subject01/Simulation');%GaussianSources');%Extents');%SNRs'); 
% [sStudy, iStudy] =bst_get('StudyWithCondition','Subject01/mind004_050924_median01_raw_clean_notch');
index = 1;
bst_call(@export_matlab, {char(sStudy.Data(1,index).FileName)},'data');
%=======================Import the LFM and Cortex==================================%
[sSurface, iSurface] = bst_get('SurfaceFileByType',[],'Cortex');
% [sSurface, iSurface] = bst_get('SurfaceFileByType',2,'Cortex');
bst_call(@export_matlab, {char(sSurface.FileName)},'Cortex');
Atlas = Cortex.Atlas(2);

[sHeadModel] = bst_get('HeadModelForStudy', iStudy);
bst_call(@export_matlab, {char(sHeadModel.FileName)},'model');
Gain=model.Gain(channelselect,:);
GridLoc=model.GridLoc;
GridOrient=model.GridOrient;

Gain = bst_gain_orient(Gain,GridOrient);
% load('MEGLeadfield.mat')
clear GridOrient
% load 'MEGLeadfield.mat';
% load 'CortexMEG.mat';
[nSensor,nSource] = size(Gain);
%% Reducing the leadfield matrix
%  L = Gain;
%  u = spm_svd(Gain*Gain');
%  L = u'*Gain;
% B = u'*B;
% clear u
%% Make output directory structure for imaging results (if doesn't already exist)
if VariousSNRs
   scenario = 'various SNRs';%'various SNRs';
   SNR1 = [-10,-5,0,5,10];
   SNIR1 = zeros(5,1)+5;
   condition = SNR1';
   K = ones(5,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousSNIRs
   scenario = 'various SNIRs';
   SNR1 = zeros(5,1)+5;
   SNIR1 = [-10,-5,0,5,10];
   condition = SNIR1';
   K = ones(5,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousExtents
   scenario = 'various extents';
   SNR1 = 5*ones(5,1);
   SNIR1 = 5*ones(5,1);
   condition = [1:5]';
   K = ones(5,1);
   DefinedArea = [2 5 10 18 32]'*1e-4*ones(1,max(K));%[0.5 4 8 14 22 32]'*1e-4*ones(1,2);% 38 48]'*1e-4;
elseif VariousChannels
   scenario = 'various channels';
   SNR1 = 5*ones(4,1);
   SNIR1 = 5*ones(4,1);
   condition = [62, 46, 32, 16]';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousPatches
   scenario = 'various patches';
   SNR1 = 5*ones(4,1);
   SNIR1 = 5*ones(4,1);
   condition = [1:4]';
   K = [1,2,3,4];
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif Test
    algorithms = {'SSSI-L2p'}%'VSSI-sL2p','VSSI-L2p'};%'VB-SCCD','SSSI-L2p','VSSI-GGD'
    ResultsLoading = [0 0 0];
    scenario = 'test';
    SNR1 = 0;
    SNIR1 = 0;
    condition = [1];
    K = 1;
    DefinedArea = 5*1e-4;
end

outpath = 'E:\项目代码(解压)\code-VSSI-Lp\VMNI_2\result_4\';
% outpath = 'E:\项目代码(解压)\code-VSSI-Lp\VMNI_2\result_4\';
for i = 1 : size(condition,1)
    path{i} = fullfile(outpath,scenario,'\',num2str(condition(i)));
    if ~exist(path{i})
        mkdir(path{i});
    end
     if ~exist([path{i} '\' 'metrics.mat'], 'file')
         metrics = [];
         save([path{i} '\' 'metrics.mat'], 'metrics')
     end
end

%% Iteration
    dim = 0;
    Miter = 30 - dim;   
    Eccentricity = sqrt(sum(GridLoc.^2,2));
    Ec = find(Eccentricity > 70*1e-3);
    %Ec = find(Eccentricity > 70);
    EstimatedArea_Mean = zeros(Miter,5);
for iter = 1:Miter    
%     ind = randperm(size(Gain,2)); 
    ind = randperm(numel(Ec));
    Thr = zeros(size(condition,1),numel(algorithms));
for iteration = 1:size(condition,1)
    fprintf('iter = %g, iteration = %g\n', iter,iteration)
    savepath = path{iteration};
    load ([path{iteration} '\' 'metrics'])
    SD = []; DLE = []; RMSE = []; AUC = []; PRE = []; REC = [];
    if any(ResultsLoading)
        SD = metrics.SD(iter,:); DLE = metrics.DLE(iter,:); RMSE = metrics.RMSE(iter,:); AUC = metrics.AUC(iter,:); PRE = metrics.PRE(iter,:); REC = metrics.REC(iter,:); nRMSE = metrics.nRMSE(iter,:); SE = metrics.SE(iter,:);
    end
%     if VariousSNIRs && (iteration==2 || iteration==3 || iteration==4)
%         continue;
%     end
%% Generate Simulated EEG Data
SNR = SNR1(iteration);
SNIR = SNIR1(iteration);
seedvox = Ec(ind(1:K(iteration))); 
% seedvox = 1436;
%seedvox = 4917;%1585;
% seedvox = 5120;
% seedvox = [4812 3321];%
%  seedvox = 3321;% for iteration effect illustration
%  seedvox = 3170;%3020;% for deep sourcces illustration
%  seedvox = 845;%5171;%661;%2360;% for extents illustration
%  seedvox = [5363 1061];%4997]%[4997 1275];% for SNR illustration
% seedvox = ind(1:K(iteration));
%        tau = [0.3 0.35 0.4 0.6]*max(Time(Activetime));%[0.25 0.29 0.4 0.50]*max(Time(Activetime));%[0.25 0.29 0.4 0.50]*max(Time(Activetime));
%        omega = [0.08 0.08 0.1 0.06];%[0.08 0.12 0.1 0.1];%*max(Time(Activetime));

%        tau = [0.1 0.15 0.3 0.4];% for imaigng purpose
%        tau = [0.25 0.29 0.4 0.50];%*max(Time(Activetime));  % for statistical simulations
tau = [0.1 0.35 0.5 0.6];omega = [0.1 0.15 0.15 0.15];%[0.07 0.05 0.10 0.10];%[0.035 0.035 0.035 0.035];%*max(Time(Activetime));
%tau = [0.1 0.2 0.5 0.6];omega = [0.1 0.15 0.15 0.15];%[0.07 0.05 0.10 0.10];%[0.035 0.035 0.035 0.035];%*max(Time(Activetime));
f = [10 11 8 9];%10*ones(1,4);%5 (Hz);
Amp = 1e-8;
StimTime = find(abs(data.Time) == min(abs(data.Time)));
TimeLen = 300;
%Time = data.Time(StimTime-0.5*TimeLen:StimTime+0.5*TimeLen-1);
Time = data.Time(1:300)
%Time = data.Time(1:40)

OPTIONS.DefinedArea    = DefinedArea(iteration,:);
OPTIONS.seedvox        = seedvox;
OPTIONS.frequency      = f;
OPTIONS.tau            = tau;
OPTIONS.omega          = omega;
OPTIONS.Amp            = Amp;
OPTIONS.GridLoc        = GridLoc;
if VariousPatches
OPTIONS.MixMatrix = [1    0    0     0;
                     0    1    0     0;
                     0    0    1     0;
                     0    0    0     1;
                     0    0    .5     .5];
elseif VariousCorrelation
    xi = condition(iteration);
%     OPTIONS.MixMatrix = [1-xi  xi   0    0;
%                          0     1    0    0;
%                          0  0   0    1;
%                          0  0   1    0];
OPTIONS.MixMatrix = [1            0            0          0;
                     xi     sqrt(1-xi^2)       0          0;
                     xi           0        sqrt(1-xi^2)   0;
                     0            0            1          0];
end
OPTIONS.uniform       = Uniform;
OPTIONS.WGN           = WGN;
OPTIONS.SNR           = SNR;
OPTIONS.SNIR          = SNIR;
OPTIONS.ar            = 0;
OPTIONS.params(:,:,1) = [ 0.8    0    0 ;
                            0  0.9  0.5 ;
                          0.4    0  0.5];

OPTIONS.params(:,:,2) = [-0.5    0    0 ;
                            0 -0.8    0 ;
                            0    0 -0.2];

OPTIONS.noisecov      = [ 0.3    0    0 ;
                            0    1    0 ;
                            0    0  0.2];
OPTIONS.SinglePoint   = 0;

if ~any(ResultsLoading)
    [Data,s_real,Result] = Simulation_Data_Generate (Gain,Cortex,Time,OPTIONS);
    ActiveVoxSeed = Result.ActiveVoxSeed;
else
    % %================= Recording from previous Results ====================%
    load ([savepath '\' 'result' num2str(iter+dim)]);
    Data = Result.B;
    s_real = Result.real;
    seedvox = Result.seedvox;
    Result.B = Data;
    [~, VertArea] = tess_area(Cortex.Vertices, Cortex.Faces);
    AreaDef = DefinedArea(iteration,:);
    ActiveVox = [];
    for k = 1:numel(seedvox)
        ActiveVoxSeed{k} = PatchGenerate(seedvox(k),Cortex.VertConn,VertArea,AreaDef(k));
        ActiveVox = union(ActiveVoxSeed{k},ActiveVox);
    end
    StimTime = size(Data,2)/2 + 1;
end
% %=======================================================================%
fprintf('Actual SNR is %g\n',20*log10(norm(Gain*s_real,'fro')/norm(Data-Gain*s_real,'fro')));
corr(s_real(seedvox,StimTime+1:end)','type','Pearson')
%% ground truth connectivity
% low_freq = 1;
% high_freq = 30;
% p = 5;
% fs = 1./(Time(2)-Time(1));
% [gamma2] = DTF(s_real(seedvox,:)',low_freq,high_freq,p,fs);
% [gamma2_sig] = DTFsigvalues(s_real(seedvox,:)', low_freq, high_freq, p, fs,[],[],[]);%, 1000, 0.05, 1);
% [gamma2] = DTFsigtest(gamma2,gamma2_sig);
% conn = zeros(size(gamma2(1)));
% for i = 1:size(gamma2,3)
%     conn = conn + gamma2(:,:,i);
% end
% conn = conn ./ size(gamma2,3);
% imagesc(conn)
% Normalize = 1;
%% Data scale
Scale = 1;
if Scale
    ScaleType = 0;  % ScaleType = 1, Simultaneously scale the MEG data and leadfiled matrix; ScaleType = 0, only scale the MEG data
    ratio = 1e-6%1e-6;%Ratio(iiter);
    if ScaleType == 0
        B = Data./ratio;
        Gain_scale = Gain;
    else
        B  = Data./ratio;
        Gain_scale = Gain./ratio;
        ratio = 1;
    end
else
    ratio = 1;
    B = Data;
    Gain_scale = Gain;
end
%% 
%======================================%
%         Leadfield Matrix normalization
%=====================================%
if LFNormlization
    LfvW = sqrt(sum(Gain_scale.^2,1).^0.3);
    Gain_scale = Gain.*kron(ones(nSensor,1),1./LfvW);
end
%% Reducing the leadfield matrix
%  u = spm_svd(Gain_scale*Gain_scale');
%  Gain_scale = u'*Gain_scale;
%  B = u'*B;
%  clear u
%% Import the Simulated EEG into Brainstorm
%     data.F = zeros(size(data.F,1),numel(Time)); 
%     data.F(channelselect,:) = B;
%     data.Time = Time;
%     %Node_Import_YZL(char(sStudy.Data(1,index).FileName),iStudy,'data');
%% Whiten measurements and lead field matrix
Bpre = B(:,1:StimTime);
Cov_n = NoiseEstimate(B,StimTime);
NoiseMethod = 'median';%'reg';%'none';%{'shrink', 'reg', 'diag', 'none', 'median'};
FourthMoment = Bpre*Bpre'./StimTime;nSamples = StimTime;
NoiseReg = 0.1;
[Cov,W] = truncate_and_regularize_covariance(Cov_n,NoiseMethod,'EEG',NoiseReg,FourthMoment,nSamples);
L = W*Gain_scale;
B = W*B;
if ~any(ResultsLoading)
    Result.Whiter = W;
end
clear W;
%%  SVD for TBFs
[Dic] = TBFSelection(B,0);%,'threshold','Permutation');'Kaiser'
% FigPath = 'E:\Figures\GLM\TBFs';
% saveas(gcf,[FigPath,'Extents.fig'])
% export_fig ([FigPath,'Extents'], '-eps', '-transparent', '-CMYK');
%% Source Estimation
% [~,MetricsInterval] = sort(mean(s_real(seedvox,StimTime:end).^2,1),'descend');
%MetricsInterval = MetricsInterval(1:10);
MetricsInterval = [];

%% Channel select
% if VariousChannels
%     if ~any(ResultsLoading)
%         cnumber = condition(iteration);
%         load(['C:\Users\guest0\Desktop\GGD-fig\channel\',num2str(cnumber),'c.mat']);
%         channels(channels>33 & channels<43) = channels(channels>33 & channels<43) - 1;
%         channels(channels>43) = channels(channels>43) - 2;
%         B = B(channels,:);
%         L = L(channels,:);
%         Result.channels = channels;
%     else
%         channels = Result.channels;
%         B = B(channels,:);
%         L = L(channels,:);
% 
%     end
% end
if VariousChannels                                                                    %这是对通道数条件的试验----通道改变的其实就是数据量
    if ~any(ResultsLoading)
        acnumber = size(B,1);                                                                    % acnumber代表的是通道数
        cnumber = condition(iteration);                                                                    %这是选择哪一种情况
        channels = sort(randsample(acnumber,cnumber));                                                                       %这是随机样本，最后进行排序
        B = B(channels,:);
        L = L(channels,:);  %L:mxn---m电极通道数  n源数   这里B和L的B = B(channels,:)因为第一维度都是通道数---m  所以第一维度取channels
        Result.channels = channels;   %Result是包含源和EEG信号灯所有的信息，是在[Data,s_real,Result] = Simulation_Data_Generate(Gain,Cortex,Time,OPTIONS);得到的
    else
        channels = Result.channels;
        B = B(channels,:);
        L = L(channels,:);
    end
end

%% VB-SCCD
Weight = logspace(-4,0,20);
% Sourcenorm = zeros(10,1);variationnorm = zeros(10,1);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;%Weight(14);%Weight(8);%0.05;%%0.15;
% opts.laplacian = 0.8;
if any(strcmpi('VB-SCCD', algorithms))
    MethodIndex = find(strcmpi('VB-SCCD', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    if ~ResultsLoading(MethodIndex)
        [S_VBSCCD] = VB_SCCD(B*Dic',L,Cortex.VertConn,opts.sparse,'transform',variation,'p',1,'tol',1e-6);
        S_vbsccd = S_VBSCCD{end}*Dic;
        if LFNormlization
            S_vbsccd = repmat(1./LfvW',1,size(S_vbsccd,2)).*S_vbsccd;
        end
        S_vbsccd = S_vbsccd*ratio;
    else
        S_vbsccd = Result.VBSCCD;
    end
%         lam_l1 = norm(B-L*s_real,'fro')/sum(abs(VariationEdge(Cortex.VertConn)*s_real))
%         [S_VBSCCD] = Lp_ADMM(B,L,Cortex.VertConn,opts.sparse,'transform',variation,'tol',1e-4,'roupar',1e17,'lam',lam_l1,'p',1,'methods',methods);
%         S_vbsccd = S_VBSCCD{end};
%         if LFNormlization
%             S_vbsccd = repmat(1./LfvW',1,size(S_vbsccd,2)).*S_vbsccd;
%         end
%         S_vbsccd = S_vbsccd*ratio;
%     else
%         S_vbsccd = Result.VBSCCD;
%     end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_vbsccd(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_vbsccd(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_vbsccd,'fro')^2/norm(B*ratio,'fro')^2;
    Result.VBSCCD = S_vbsccd;
end 

%% SISSY
Weight = logspace(-4,0,20);
% Sourcenorm = zeros(10,1);variationnorm = zeros(10,1);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;%Weight(14);%Weight(8);%0.05;%%0.15;
% opts.laplacian = 0.8;
if any(strcmpi('SISSY', algorithms))
    MethodIndex = find(strcmpi('SISSY', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    if ~ResultsLoading(MethodIndex)
        [S_SISSY] = SISSY(B*Dic',L,Cortex.VertConn,opts.sparse,'transform',variation,'p',1,'tol',1e-6);
        S_sissy = S_SISSY{end}*Dic;
        if LFNormlization
            S_sissy = repmat(1./LfvW',1,size(S_sissy,2)).*S_sissy;
        end
        S_sissy = S_sissy*ratio;
    else
        S_sissy = Result.SISSY;
    end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_sissy(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_sissy(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_sissy,'fro')^2/norm(B*ratio,'fro')^2;
    Result.SISSY = S_sissy;
end


%% VSSI-L1N
Weight = logspace(-4,0,20);
% Sourcenorm = zeros(10,1);variationnorm = zeros(10,1);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;%Weight(14);%Weight(8);%0.05;%%0.15;
% opts.laplacian = 0.8;
if any(strcmpi('VSSI-L1N', algorithms))
    MethodIndex = find(strcmpi('VSSI-L1N', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    methods = algorithms{MethodIndex};

    if ~ResultsLoading(MethodIndex)
        p = 0.8;
%         Mtemp = [opts.sparse*sparse(1:size(s_real,1),1:size(s_real,1),1);VariationEdge(Cortex.VertConn)];
%         lam = norm(B./ratio_Lp-L*s_real,'fro')/sum(sqrt(sum((Mtemp*s_real).^2,2)).^p)
        [S_1N] = VSSI_L1N(B*Dic',L,Cortex.VertConn,opts.sparse,'transform',variation,'p',p,'tol',1e-6,'Dic',Dic);
        S_vssil1n = S_1N{end}*Dic;
        if LFNormlization
            S_vssil1n = repmat(1./LfvW',1,size(S_vssil1n,2)).*S_vssil1n;
        end
        S_vssil1n = S_vssil1n*ratio;
    else
        S_vssil1n = Result.VSSIL1N;
    end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_vssil1n(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_vmnil21n(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_vssil1n,'fro')^2/norm(B*ratio,'fro')^2;
    Result.VSSIL1N = S_vssil1n;
end 
  
  %% L_2 norm inverse solver
  if any(strcmpi('wMNE', algorithms))
    MethodIndex = find(strcmpi('wMNE', algorithms)~=0);
    if ~ResultsLoading(MethodIndex)
        TwMNE = MNE(B,model.Gain(channelselect,:),L,Cortex,'wMNE');
        if LFNormlization
            TwMNE = repmat(1./LfvW',1,size(TwMNE,2)).*TwMNE;
        end
        S_wMNE = TwMNE*B*ratio;
    else
        TwMNE = Result.wMNE;
        S_wMNE = TwMNE*Result.Whiter*Result.B;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_wMNE(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_wMNE(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_wMNE,'fro')^2/norm(B*ratio,'fro')^2;
    Result.wMNE = TwMNE; 
  end

 
 if any(strcmpi('LORETA', algorithms))
    MethodIndex = find(strcmpi('LORETA', algorithms)~=0);
    if ~ResultsLoading(MethodIndex)
        TLORETA = MNE(B,model.Gain(channelselect,:),L,Cortex,'LORETA');
        if LFNormlization
            TLORETA = repmat(1./LfvW',1,size(TLORETA,2)).*TLORETA;
        end
        S_LORETA = TLORETA*B*ratio;
    else
        TLORETA = Result.LORETA;
        S_LORETA = TLORETA*Result.Whiter*Result.B;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_LORETA(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_LORETA(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_LORETA,'fro')^2/norm(B*ratio,'fro')^2;
    Result.LORETA = TLORETA;
 end
  
 if any(strcmpi('sLORETA', algorithms))
    MethodIndex = find(strcmpi('sLORETA', algorithms)~=0);
    if ~ResultsLoading(MethodIndex)
        TsLORETA = MNE(B,model.Gain(channelselect,:),L,Cortex,'sloreta');
        if LFNormlization
            TsLORETA = repmat(1./LfvW',1,size(TsLORETA,2)).*TsLORETA;
        end
        S_sLORETA = TsLORETA*B*ratio;
    else
        TsLORETA = Result.sLORETA;
        S_sLORETA = TsLORETA*Result.Whiter*Result.B;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_sLORETA(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_sLORETA(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    Result.sLORETA = TsLORETA;
 end
 
  if any(strcmpi('dspm', algorithms))
    MethodIndex = find(strcmpi('dspm', algorithms)~=0);
    if ~ResultsLoading(MethodIndex)
        Tdspm = MNE(B,model.Gain(channelselect,:),L,Cortex,'dspm','reg',0);
        if LFNormlization
            Tdspm = repmat(1./LfvW',1,size(Tdspm,2)).*Tdspm;
        end
        S_dspm = Tdspm*B*ratio;
    else
        Tdspm = Result.dspm;
        S_dspm = Tdspm*Result.Whiter*Result.B;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_dspm(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_dspm(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    Result.dspm = Tdspm;
  end
 
  if any(strcmpi('MNE', algorithms))
    MethodIndex = find(strcmpi('MNE', algorithms)~=0);
    if ~ResultsLoading(MethodIndex)
        TMNE = MNE(B,model.Gain(channelselect,:),L,Cortex,'mne');
        if LFNormlization
            TMNE = repmat(1./LfvW',1,size(TMNE,2)).*TMNE;
        end
        S_MNE = TMNE*B*ratio;
    else
        TMNE = Result.MNE;
        S_MNE = TMNE*Result.Whiter*Result.B;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_MNE,s_real,ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_MNE(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_MNE,'fro')^2/norm(B*ratio,'fro')^2;
    Result.MNE = TMNE;
  end
 %% SBL
 if any(strcmpi('Champagne', algorithms))
     MethodIndex = find(strcmpi('Champagne', algorithms)~=0);
     if ~ResultsLoading(MethodIndex)
         [SBLkernel,par_SBL]=SBL(B,L,'epsilon',1e-4,'flags',1,'prune',[1,1e-6]);
         S_SBL = SBLkernel*B*ratio;
         if LFNormlization
             S_SBL = repmat(1./LfvW',1,size(S_SBL,2)).*S_SBL;
         end
     else
         SBLkernel = Result.SBLkernel;
         S_SBL = SBLkernel*Result.Whiter*Result.B;
     end
     
     [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
         = PerformanceMetric(GridLoc,S_SBL(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
     A = ROCextent(s_real(:,StimTime:end),S_SBL(:,StimTime:end),Cortex,seedvox,'flag',2);
     AUC(1,MethodIndex) = median(A.mean);
     EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_SBL,'fro')^2/norm(B*ratio,'fro')^2;
     Result.SBLkernel = SBLkernel;
 end
%% FISTA parameter(when regualar parameter isn't learned from CV, this part is needed)
if any(strcmpi('SSI', algorithms))
    MethodIndex = find(strcmpi('SSI', algorithms)~=0);
    addpath('build');
    setenv('MKL_NUM_THREADS','1')
    setenv('MKL_SERIAL','YES')
    setenv('MKL_DYNAMIC','NO')
    gamma_min = 8;
    gamma_max = 12;
    n_gamma = 100;
    gamma_all = logspace(gamma_min,gamma_max,n_gamma);
    
    format compact;
    param.numThreads=-1; % all cores (-1 by default)
    param.verbose=true;   % verbosity, false by default
    param.it0=50;      % frequency for duality gap computations
    param.max_it = 2000; % maximum number of iterations
    param.L0 = 0.1;
    param.tol = 1e-3;
    param.intercept  = false;
    param.pos = false;
    param.compute_gram = false;
    
    param.loss = 'square';
    fprintf('\nFISTA + Regression l1l2 \n');
    param.regul = 'l1l2';
    
    Phi = (orth(Dic'))';
    W0 = zeros(nSource,size(Phi,1));
    % param.lambda = 1e11*ratio;
    
    param.lambda = ratio*size(Phi,1)*nSensor/norm(s_real*Phi','fro');
    
    [W, optim_info] = mexFistaFlat(B*Phi',L,W0,param);
    S_sparse = W*Phi*ratio;
    if LFNormlization
        S_sparse = repmat(1./LfvW',1,size(S_sparse,2)).*S_sparse;
    end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_sparse(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_sparse(:,StimTime:end),Cortex,seedvox,'flag',2);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_sparse,'fro')^2/norm(B*ratio,'fro')^2;
    Result.Sparse = S_sparse;
end
%% BESTIES
if any(strcmpi('BESTIES', algorithms))
    MethodIndex = find(strcmpi('BESTIES', algorithms)~=0);
    [S_BESTIES,par] = BESTIES(B,L,Dic,Cortex.VertConn,'epsilon',5*1e-5,'vb_glm',0);
    if LFNormlization
        S_BESTIES = repmat(1./LfvW',1,size(S_BESTIES,2)).*S_BESTIES;
    end
    S_BESTIES = S_BESTIES*ratio;
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_BESTIES(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_BESTIES(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_BESTIES,'fro')^2/norm(B*ratio,'fro')^2;
    Result.BESTIES = S_BESTIES;
end
%% GLM
if any(strcmpi('GLM', algorithms))
    MethodIndex = find(strcmpi('GLM', algorithms)~=0);
    [S_vbglm,par_vbglm] = BESTIES(B,L,Dic,Cortex.VertConn,'vb_glm',1);
    if LFNormlization
        S_vbglm = repmat(1./LfvW',1,size(S_vbglm,2)).*S_vbglm;
    end
    S_vbglm = S_vbglm*ratio;
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_vbglm(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_vbglm(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_vbglm,'fro')^2/norm(B*ratio,'fro')^2;
    Result.vbglm = S_vbglm;
end
%% MFOCUSS
if any(strcmpi('MFOCUSS', algorithms))
    MethodIndex = find(strcmpi('MFOCUSS', algorithms)~=0);
%     lambda = ratio*size(B,2)*nSensor/sum(sum(s_real.^2,2).^(0.8/2));
    Lambda = logspace(3,8,10);
%     P = sqrt(sum((L'*B).^2,2));
%     Lambda = logspace(-8,0,20)*max(P);
% % ===================== Select the optimal regularizer =================%  

% % ======================================================================%
if ~ResultsLoading(MethodIndex)
    lambda = 1e-6*size(B,2)*nSensor/sum(sum(s_real.^2,2).^(0.8/2));%*0.05;
%     if SNR == 5 || SNR == 10 || SNR == 15
%         lambda = Lambda(5);
%     end
%     lambda = Lambda(5);%Lambda(i); %% for single patch, lambda = Lambda(7); for extent = 3, lambda = Lambda(6); SNR = 10, lambda = Lambda(5);
    [S_FOCUSS,gamma_ind,gamma_est,count] = MFOCUSS(L,B, lambda,'PRINT',1);

%     [S_FOCUSS,gamma_mean] = M_FOCUSS(L,B,2,16);

    if LFNormlization
        S_FOCUSS = repmat(1./LfvW',1,size(S_FOCUSS,2)).*S_FOCUSS;
    end
    S_FOCUSS = S_FOCUSS*ratio;
else
    S_FOCUSS = Result.FOCUSS;
end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_FOCUSS(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_FOCUSS(:,StimTime:end),Cortex,seedvox,'flag',2);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_FOCUSS,'fro')^2/norm(B*ratio,'fro')^2;
    Result.FOCUSS = S_FOCUSS;
end

%% test Beamformer
if any(strcmpi('Beamformer', algorithms))
    MethodIndex = find(strcmpi('Beamformer', algorithms)~=0);
    beamformer = 2;
    if ~ResultsLoading(MethodIndex)
        switch beamformer
            case 1
                Options.LFVWeight=0;
                [S_beamformer,Kernel,wNorm,sPower] = SourceImagingBF_BasicMethod(L,B,[],0.5,[],Options);
                
            case 2
                Options.firstEye=0;
                Options.MaxIter=100;
                Options.stopCriteria=0.05;
                Options.LFVWeight=1;
                [S_beamformer,Kernel,wNorm,sPower] = SourceImagingBF_YZLDirectIterate(L,B,[],0.2,[],Options);
                
            case 3
                Options.SNRprune=6;
                Options.firstEye=0;
                Options.MaxIter=100;
                Options.stopCriteria=0.05;
                Options.LFVWeight=1;
                [S_beamformer,Kernel,wNorm,sPower] = SourceImagingBF_YZLPrunedIterate(L,B,[],0.2,[],Options);
                
            case 4
                Options.rankRatio=0.95;
                Options.SNRprune=6;
                Options.firstEye=0;
                Options.MaxIter=100;
                Options.stopCriteria=0.05;
                Options.LFVWeight=1;
                [S_beamformer,Kernel,wNorm,sPower] = SourceImagingBF_YZLProjectedPrunedIterate(L,B,[],0.1,[],Options);
        end
        if LFNormlization
            S_beamformer = repmat(1./LfvW',1,size(S_beamformer,2)).*S_beamformer;
        end
        S_beamformer = S_beamformer*ratio;
    else
        S_beamformer = Result.Beamformer;
    end
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_beamformer(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_beamformer(:,StimTime:end),Cortex,seedvox,'flag',1);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_beamformer,'fro')^2/norm(B*ratio,'fro')^2;
    Result.Beamformer = S_beamformer;
end


%% Save Results on Disk
method = 1:numel(algorithms);
metrics.AUC(dim+iter,method) = AUC;%(method);
metrics.SD(dim+iter,method) = SD;%(method); 
metrics.DLE(dim+iter,method) = DLE;%(method);
metrics.RMSE(dim+iter,method) = RMSE;%(method);
metrics.nRMSE(dim+iter,method) = nRMSE;%(method);
metrics.SE(dim+iter,method) = SE;%(method);
metrics.PRE(dim+iter,method) = PRE;
metrics.REC(dim+iter,method) = REC;


save_file_name=[savepath '\' 'result' num2str(dim+iter)];
save (save_file_name,'Result');
save([savepath '\' 'metrics.mat'], 'metrics')

end
end
toc
save([savepath '\' 'GridLoc.mat'], 'GridLoc')
save([savepath '\' 'Cortex.mat'], 'Cortex')
save([savepath '\' 'methods.mat'], 'algorithms')
runningtime = toc;