clear
clc;
warning off

path = '';
addpath(genpath(path));

%dataName = 'ORL'; %%% flower17; proteinFold;  ucsd-mit_caltech-101-mkl
%load([path,dataName,'_Kmatrix'],'KH','Y','P');

alpha = 0.9;
numclass =5;
folder_name = "Experiment1-label=1-smoothing_alpha="+num2str(alpha)+"-norm=True";
%folderName = "Experiment_Sp-label=1-norm=True";
if ~exist(folder_name+'/rnaseq-kms', 'dir')
   mkdir(folder_name+'/rnaseq-kms');
end
if ~exist(folder_name+'/rppa-kms', 'dir')
   mkdir(folder_name+'/rppa-kms');
end
if ~exist(folder_name+'/som-kms', 'dir')
   mkdir(folder_name+'/som-kms');
end
folderNameOut = folder_name+"/labels";
if ~exist(folderNameOut, 'dir')
   mkdir(folderNameOut);
end
KH_rs_fn = unzip(folder_name + '/rnaseq-kms.npz',folder_name+'/rnaseq-kms');
KH_rp_fn = unzip(folder_name + '/rppa-kms.npz',folder_name+'/rppa-kms');
KH_som_fn = unzip(folder_name + '/som-kms.npz',folder_name+'/som-kms');

KH_rs = readNPY(folder_name + '/rnaseq-kms/kms.npy');
numker_rs = size(KH_rs,1);
KH_rp = readNPY(folder_name + '/rppa-kms/kms.npy');
numker_rp = size(KH_rp,1);
KH_som = readNPY(folder_name + '/som-kms/kms.npy');
numker_som = size(KH_som,1);

numsample = size(KH_rs,3);

KH = reshape(KH_rs(1,:,:),numsample,numsample);
for i = 2:numker_rs
    data = reshape(KH_rs(i,:,:),numsample,numsample);
    KH = cat(3,KH,data);
end

for i = 1:numker_rp
    data = reshape(KH_rp(i,:,:),numsample,numsample);
    KH = cat(3,KH,data);
end

for i = 1:numker_som
    data = reshape(KH_som(i,:,:),numsample,numsample);
    KH = cat(3,KH,data);
end

KH = kcenter(KH);
KH = knorm(KH);
%numclass = length(unique(Y));
numker = size(KH,3);
%%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma0 = ones(numker,1)/numker;
avgKer = mycombFun(KH,gamma0);
[H_normalized1] = mykernelkmeans(avgKer,numclass);
indx = untitled(H_normalized1,numclass);
name = folderNameOut+"/smspk-all-kmeans-"+int2str(numclass)+"lab";
csvwrite(name,indx)

%%%%%%%%%%---AAAI-16----%%%%%%%%
lambdaset2 = 2.^[-15:3:15];
%lambdaset2 = 2.^[5:2:6];

M = calculateM(KH);
for il =1:length(lambdaset2)
    [H_normalized2,gamma2,obj2] = myregmultikernelclustering(KH,M,numclass,lambdaset2(il));
    indx2 = untitled(H_normalized2,numclass);
    name = strcat(folderNameOut+"/smspk-all-mkkm-"+int2str(numclass)+"lab-lambda"+ num2str(il));
    csvwrite(name,indx2)
end
%res(:,2) = [max(accval2);max(nmival2);max(purval2)];

%% Contact 1022xinwang.liu@gmail.com
%%%%%%%-----Citation----------------------%%%%%%%
%% Xinwang Liu, Yong Dou, Jianping Yin, Lei Wang, En Zhu: 
%% Multiple Kernel k-Means Clustering with Matrix-Induced Regularization. AAAI 2016: 1888-1894

