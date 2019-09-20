function e = smspk_clustering_drop_fnc(folder_name,numclass, drop_percent)

warning off
numclass = double(numclass);

disp(folder_name);
folder_name = strcat(folder_name,"");
disp(~exist(strcat(folder_name+"/rnaseq-kms"), "dir"));
if ~exist(folder_name+"/rnaseq-kms", "dir")
   mkdir(folder_name+"/rnaseq-kms");
   KH_rs_fn = unzip(folder_name + "/rnaseq-kms.npz",folder_name+"/rnaseq-kms");
end
if ~exist(folder_name+"/rppa-kms", "dir")
   mkdir(folder_name+"/rppa-kms");
   KH_rp_fn = unzip(folder_name + "/rppa-kms.npz",folder_name+"/rppa-kms");
end
if ~exist(folder_name+"/som-kms", "dir")
   mkdir(folder_name+"/som-kms");
   KH_som_fn = unzip(folder_name + "/som-kms.npz",folder_name+"/som-kms");
end
folderNameOut = folder_name+"/labels_dropped"+num2str(drop_percent);
if ~exist(folderNameOut, "dir")
   mkdir(folderNameOut);
end

KH_rs = readNPY(folder_name + "/rnaseq-kms/kms.npy");
numker_rs = size(KH_rs,1);
KH_rp = readNPY(folder_name + "/rppa-kms/kms.npy");
numker_rp = size(KH_rp,1);
KH_som = readNPY(folder_name + "/som-kms/kms.npy");
numker_som = size(KH_som,1);

numsample = size(KH_rs,3);

KH_rs = permute(KH_rs, [2 3 1]);
KH_rp = permute(KH_rp, [2 3 1]);
KH_som = permute(KH_som, [2 3 1]);

KH = cat(3, KH_rs,KH_rp,KH_som);
dropped = [];
stayed = [];
total = numsample*numsample;
limit = (drop_percent*total)/100;
for k = 1:size(KH,3)
   if nnz(KH(:,:,k)) < limit
       dropped = [dropped k];
   else
       stayed = [stayed k];
   end
end
KH(:,:,dropped) = [];
stayed = stayed';
name_dropped = strcat(folderNameOut+"/dropped_kernels");
csvwrite(name_dropped, dropped)

KH = kcenter(KH);
KH = knorm(KH);
%numclass = length(unique(Y));
numker = size(KH,3);
gamma0 = ones(numker,1)/numker;
avgKer = mycombFun(KH,gamma0);
M = calculateM(KH);
%for numclass = [2, 3, 4, 5]
%%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H_normalized1] = mykernelkmeans(avgKer,numclass);
indx = untitled(H_normalized1,numclass);
name = folderNameOut+"/smspk-all-kmeans-"+int2str(numclass)+"lab";
csvwrite(name,indx)

%%%%%%%%%%---AAAI-16----%%%%%%%%
lambdawriter = [-15:3:15];
lambdaset2 = 2.^[-15:3:15];
%lambdaset2 = 2.^[5:2:6];
%lambdaset2 = 2.^[0];
for il =1:length(lambdaset2)
    [H_normalized2,gamma2,obj2] = myregmultikernelclustering(KH,M,numclass,lambdaset2(il));
    indx2 = untitled(H_normalized2,numclass);
    name = strcat(folderNameOut+"/smspk-all-mkkm-"+int2str(numclass)+"lab-log(lambda)="+ num2str(lambdawriter(il)));
    csvwrite(name,indx2)
    name_gamma = strcat(folderNameOut+"/smspk-all-mkkm-"+int2str(numclass)+"lab-log(lambda)="+ num2str(lambdawriter(il))+"weights");
    gammas = [stayed gamma2];
    csvwrite(name_gamma,gammas)
end
e=drop_percent;
end

%% Contact 1022xinwang.liu@gmail.com
%%%%%%%-----Citation----------------------%%%%%%%
%% Xinwang Liu, Yong Dou, Jianping Yin, Lei Wang, En Zhu: 
%% Multiple Kernel k-Means Clustering with Matrix-Induced Regularization. AAAI 2016: 1888-1894

