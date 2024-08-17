clc
clear
tic
%%
algorithms = {'VSSI-L1N','SISSY','VB-SCCD','wMNE','LORETA'};
ResultsLoading = [0 0 0 0 0];
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
else
    channelselect=[1:32 34:42 44:59 61:63]; % for Real noise simulations
end
[sStudy, iStudy] = bst_get('StudyWithCondition','Subject01/Simulation');
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

outpath = '';
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
    EstimatedArea_Mean = zeros(Miter,5);
for iter = 1:Miter    
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
tau = [0.1 0.35 0.5 0.6];omega = [0.1 0.15 0.15 0.15];
Amp = 1e-8;
StimTime = find(abs(data.Time) == min(abs(data.Time)));
TimeLen = 300;
Time = data.Time(1:300)

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
%% Data scale
Scale = 1;
if Scale
    ScaleType = 0;
    ratio = 1e-6;
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
%% Whiten measurements and lead field matrix
Bpre = B(:,1:StimTime);
Cov_n = NoiseEstimate(B,StimTime);
NoiseMethod = 'median';
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
if VariousChannels                                                                 
    if ~any(ResultsLoading)
        acnumber = size(B,1);
        cnumber = condition(iteration);
        channels = sort(randsample(acnumber,cnumber));
        B = B(channels,:);
        L = L(channels,:);
        Result.channels = channels;
    else
        channels = Result.channels;
        B = B(channels,:);
        L = L(channels,:);
    end
end

%% VB-SCCD
Weight = logspace(-4,0,20);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;
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
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_vbsccd(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_vbsccd(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_vbsccd,'fro')^2/norm(B*ratio,'fro')^2;
    Result.VBSCCD = S_vbsccd;
end 

%% SISSY
Weight = logspace(-4,0,20);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;
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
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;
if any(strcmpi('VSSI-L1N', algorithms))
    MethodIndex = find(strcmpi('VSSI-L1N', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    methods = algorithms{MethodIndex};

    if ~ResultsLoading(MethodIndex)
        p = 0.8;
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
 
%% Save Results on Disk
method = 1:numel(algorithms);
metrics.AUC(dim+iter,method) = AUC;
metrics.SD(dim+iter,method) = SD;
metrics.DLE(dim+iter,method) = DLE;
metrics.RMSE(dim+iter,method) = RMSE;
metrics.nRMSE(dim+iter,method) = nRMSE;
metrics.SE(dim+iter,method) = SE;
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