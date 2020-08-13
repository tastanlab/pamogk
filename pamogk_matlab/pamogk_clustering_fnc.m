function e = pamogk_clustering_fnc(folder_name, numclass, drop_percent)

warning off
numclass = double(numclass);

ROOT_DIR = strcat(folder_name,"");

if ~exist(ROOT_DIR+"/rnaseq-kms", "dir")
   mkdir(ROOT_DIR+"/rnaseq-kms");
   KH_rs_fn = unzip(ROOT_DIR + "/rnaseq-kms.npz",ROOT_DIR+"/rnaseq-kms");
end
if ~exist(ROOT_DIR+"/rppa-kms", "dir")
   mkdir(ROOT_DIR+"/rppa-kms");
   KH_rp_fn = unzip(ROOT_DIR + "/rppa-kms.npz",ROOT_DIR+"/rppa-kms");
end
if ~exist(ROOT_DIR+"/som-kms", "dir")
  mkdir(ROOT_DIR+"/som-kms");
  KH_som_fn = unzip(ROOT_DIR + "/som-kms.npz", ROOT_DIR+"/som-kms");
end

folderNameOut = ROOT_DIR+"/labels_dropped_"+num2str(drop_percent);
if ~exist(folderNameOut, "dir")
   mkdir(folderNameOut); 
end

KH_rs = readNPY(ROOT_DIR + "/rnaseq-kms/kms.npy");
numker_rs = size(KH_rs,1);
KH_rp = readNPY(ROOT_DIR + "/rppa-kms/kms.npy");
numker_rp = size(KH_rp,1);
KH_som = readNPY(ROOT_DIR + "/som-kms/kms.npy");
numker_som = size(KH_som,1);

numsample = size(KH_rs,3);

%KH_rs = permute(KH_rs, [2 3 1]);
%KH_rp = permute(KH_rp, [2 3 1]);
%KH_som = permute(KH_som, [2 3 1]);

KH = cat(1, KH_rs,KH_rp,KH_som);
%KH = KH_rp;
dropped = [];
stayed = [];
KH_SNF_cell = {};
total = numsample*numsample;
limit = (drop_percent*total)/100;
for k = 1:size(KH,1)
   if nnz(KH(k,:,:)) < limit
       dropped = [dropped k];
   else
       stayed = [stayed k];
       new_kernel = reshape(KH(k,:,:),size(KH,2), size(KH,3));
       new_kernel = kcenter(new_kernel);
       new_kernel = knorm(new_kernel);
       KH_SNF_cell = [KH_SNF_cell,new_kernel];
   end
end
KH(dropped,:,:) = [];
stayed = stayed';
name_dropped = strcat(folderNameOut+"/dropped_kernels");
csvwrite(name_dropped, dropped)

K = 20;%number of neighbors, usually (10~30)
T = 20; %Number of Iterations, usually (10~20)

snf_kmeans_loc = folderNameOut+"/pamogk-snf-kmeans-"+int2str(numclass)+"lab";
snf_spectral_loc = folderNameOut+"/pamogk-snf-spectral-"+int2str(numclass)+"lab";

W = SNF(KH_SNF_cell, K, T);
[H_normalized_snf] = mykernelkmeans(W, numclass);
km_snf = normalized_kmeans(H_normalized_snf, numclass);
csvwrite(snf_kmeans_loc,km_snf)

group = SpectralClustering(W, numclass);
csvwrite(snf_spectral_loc,group)

kmeans_name = folderNameOut+"/pamogk-kmeans-"+int2str(numclass)+"lab";

KH= permute(KH, [2 3 1]);
KH = kcenter(KH);
KH = knorm(KH);
%numclass = length(unique(Y));
numker = size(KH,3);
gamma0 = ones(numker,1)/numker;
avgKer = mycombFun(KH,gamma0);
M = calculateM(KH);
%for numclass = [2, 3, 4, 5]
%%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kmeans_name = folderNameOut+"/pamogk-kmeans-"+int2str(numclass)+"lab";
[H_normalized1] = mykernelkmeans(avgKer,numclass);
indx = normalized_kmeans(H_normalized1,numclass);
csvwrite(kmeans_name,indx)

%%%%%%%%%%---AAAI-16----%%%%%%%%
lambdaset2 = 2.^(-15:3:15);
lambdasetname = -15:3:15;
%lambdaset2 = 2.^[5:2:6];
%lambdaset2 = 2.^[0];
for il =1:length(lambdaset2)
    [H_normalized2,gamma2,obj2] = myregmultikernelclustering(KH,M,numclass,lambdaset2(il));
    indx2 = normalized_kmeans(H_normalized2,numclass);
    name = strcat(folderNameOut+"/pamogk-mkkm-"+int2str(numclass)+"lab-loglambda="+ num2str(lambdasetname(il)));
    csvwrite(name,indx2)
    name_gamma = strcat(folderNameOut+"/pamogk-mkkm-"+int2str(numclass)+"lab-loglambda="+ num2str(lambdasetname(il))+"weights");
    gammas = [stayed gamma2];
    csvwrite(name_gamma,gammas)
    name_obj = strcat(folderNameOut+"/pamogk-mkkm-"+int2str(numclass)+"lab-loglambda="+ num2str(lambdasetname(il))+"obj");
    csvwrite(name_obj,obj2)
end
%end
e=drop_percent;
end

%res(:,2) = [max(accval2);max(nmival2);max(purval2)];

%% Contact 1022xinwang.liu@gmail.com
%%%%%%%-----Citation----------------------%%%%%%%
%% Xinwang Liu, Yong Dou, Jianping Yin, Lei Wang, En Zhu: 
%% Multiple Kernel k-Means Clustering with Matrix-Induced Regularization. AAAI 2016: 1888-1894

