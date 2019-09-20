function e = stmk_clustering_drop_fnc(folder_name,numclass, drop_percent)

%clear
%clc;
warning off
numclass = double(numclass);

disp(folder_name);
folder_name = strcat(folder_name,"");
disp(~exist(strcat(folder_name+"/rs-o-kms"), "dir"));
if ~exist(folder_name+'/rs-o-kms', 'dir')
   mkdir(folder_name+'/rs-o-kms');
   KH_rs_fn = unzip(folder_name + '/rs-o-kms.npz',folder_name+'/rs-o-kms');
end
if ~exist(folder_name+'/rs-u-kms', 'dir')
   mkdir(folder_name+'/rs-u-kms');
   KH_rs_fn = unzip(folder_name + '/rs-u-kms.npz',folder_name+'/rs-u-kms');
end
if ~exist(folder_name+'/rp-o-kms', 'dir')
   mkdir(folder_name+'/rp-o-kms');
   KH_rs_fn = unzip(folder_name + '/rp-o-kms.npz',folder_name+'/rp-o-kms');
end
if ~exist(folder_name+'/rp-u-kms', 'dir')
   mkdir(folder_name+'/rp-u-kms');
   KH_rs_fn = unzip(folder_name + '/rp-u-kms.npz',folder_name+'/rp-u-kms');
end
if ~exist(folder_name+'/som-kms', 'dir')
   mkdir(folder_name+'/som-kms');
   KH_som_fn = unzip(folder_name + '/som-kms.npz',folder_name+'/som-kms');
end
folderNameOut = folder_name+"/labels_dropped"+num2str(drop_percent);
if ~exist(folderNameOut, 'dir')
   mkdir(folderNameOut);
end
name = folderNameOut+"/stmk-all-kmeans-"+int2str(numclass)+"lab";
if exist(name, 'file')
   e=0
   return
end

KH_rs_o = readNPY(folder_name + '/rs-o-kms/kms.npy');
numker_rs_o = size(KH_rs_o,1);
KH_rs_u = readNPY(folder_name + '/rs-u-kms/kms.npy');
numker_rs_u = size(KH_rs_u,1);
KH_rp_o = readNPY(folder_name + '/rp-o-kms/kms.npy');
numker_rs_o = size(KH_rp_o,1);
KH_rp_u = readNPY(folder_name + '/rp-u-kms/kms.npy');
numker_rs_u = size(KH_rp_u,1);
KH_som = readNPY(folder_name + '/som-kms/kms.npy');
numker_som = size(KH_som,1);

numsample = size(KH_rs_o,3);

KH_rs_o = permute(KH_rs_o, [2 3 1]);
KH_rp_o = permute(KH_rp_o, [2 3 1]);
KH_rs_u = permute(KH_rs_u, [2 3 1]);
KH_rp_u = permute(KH_rp_u, [2 3 1]);
KH_som = permute(KH_som, [2 3 1]);

KH = cat(3, KH_rs_o,KH_rs_u,KH_rp_o,KH_rp_u,KH_som);
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
%%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H_normalized1] = mykernelkmeans(avgKer,numclass);
indx = untitled(H_normalized1,numclass);
name = folderNameOut+"/stmk-all-kmeans-"+int2str(numclass)+"lab";
csvwrite(name,indx)

%%%%%%%%%%---AAAI-16----%%%%%%%%
lambdaset2 = 2.^[-15:3:15];
%lambdaset2 = 2.^[5:2:6];
%lambdaset2 = 2.^[0];
for il =1:length(lambdaset2)
    [H_normalized2,gamma2,obj2] = myregmultikernelclustering(KH,M,numclass,lambdaset2(il));
    indx2 = untitled(H_normalized2,numclass);
    name = strcat(folderNameOut+"/stmk-all-mkkm-"+int2str(numclass)+"lab-lambda"+ num2str(il));
    csvwrite(name,indx2)
    name_gamma = strcat(folderNameOut+"/stmk-all-mkkm-"+int2str(numclass)+"lab-lambda"+ num2str(il)+"weights");
    gammas = [stayed gamma2];
    csvwrite(name_gamma,gammas)
end
%end
e=drop_percent;
end

%res(:,2) = [max(accval2);max(nmival2);max(purval2)];

%% Contact 1022xinwang.liu@gmail.com
%%%%%%%-----Citation----------------------%%%%%%%
%% Xinwang Liu, Yong Dou, Jianping Yin, Lei Wang, En Zhu: 
%% Multiple Kernel k-Means Clustering with Matrix-Induced Regularization. AAAI 2016: 1888-1894

