clear all;
warning off;
clc;
bits = 128;
run = 3;
model.MAP=[];

%% Dataset Loading
% load mir_cnn.mat;
load('I:\zcq\mir_cnn.mat')
fprintf('MIR Flickr_CNN dataset loaded...\n');
%% Parameters Setting
alpha = 1e-3;
beta = 1e-1;
lambda = 1e-1;
eta = 1e-1;
SL = 3;
SU = 1;

%% centralization
fprintf('centralizing data...\n');
Ntrain = size(I_tr,1);
% get anchors
n_anchors = 300;
sample = randsample(Ntrain, n_anchors);
anchorI = I_tr(sample,:);
anchorT = T_tr(sample,:);
sigmaI=100;
sigmaT=100;%MIF FLichr and NUS WIDE
% sigmaI=36;
% sigmaT=4;% MS COCO
Phi_trainI = exp(-sqdist(I_tr,anchorI)/(2*sigmaI*sigmaI));
Phi_trainI = [Phi_trainI, ones(size(Phi_trainI,1),1)];
Pht_trainT = exp(-sqdist(T_tr,anchorT)/(2*sigmaT*sigmaT));
Pht_trainT = [Pht_trainT, ones(size(Pht_trainT,1),1)];

Nlab = 5000; % the number of the labeled data
L = L_tr(1:Nlab,:);
S = L*L';
S(S>=1) = bits;
S(S==0) = -bits;
Y = L;

I_temp = Phi_trainI;
T_temp = Pht_trainT;
param.SU = SU;
param.SL = SL;
param.alpha = alpha;
param.beta = beta;
param.lambda = lambda ;
param.eta = eta;
param.bits = bits;
%% Training & Evaluation Process
fprintf('\n============================================Start training ESMH============================================\n');
for j = 1:run
    % Training model
    [F,W1, W2,B,G,mu1, mu2] = solve_ESMH(I_temp, T_temp, S,Y, param);
    
    % Evaluation
    fprintf('Evaluation...\n');
    Phi_dbI = exp(-sqdist(I_db,anchorI)/(2*sigmaI*sigmaI));
    Phi_dbI = [Phi_dbI, ones(size(Phi_dbI,1),1)];
    Pht_dbT = exp(-sqdist(T_db,anchorT)/(2*sigmaT*sigmaT));
    Pht_dbT = [Pht_dbT, ones(size(Pht_dbT,1),1)];
    Phi_testI = exp(-sqdist(I_te,anchorI)/(2*sigmaI*sigmaI));
    Phi_testI = [Phi_testI, ones(size(Phi_testI,1),1)];
    Pht_testT = exp(-sqdist(T_te,anchorT)/(2*sigmaT*sigmaT));
    Pht_testT = [Pht_testT, ones(size(Pht_testT,1),1)];
    
    B_db = mu1*Phi_dbI*W1+mu2*Pht_dbT*W2>0;
    B_test = mu1*Phi_testI*W1+mu2*Pht_testT*W2>0;
    Vdb = compactbit(B_db);
    Vtest = compactbit(B_test);
    
    Dhamm = hammingDist(Vdb, Vtest);
    [MAP] = perf_metric4Label(L_db, L_te, Dhamm);
    map(j) = MAP;
    fprintf('============================================The number of the labeled data is %d =============================================\n', Nlab);
    fprintf('============================================%d bits ESMH mAP over %d iterations:%.4f=============================================\n', param.bits, run, mean(map));
end
