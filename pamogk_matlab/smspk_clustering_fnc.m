function e = pamogk_clustering_fnc(kernel_path, result_dir, numclass, drop_percent, typ)

warning off
numclass = double(numclass);

load(kernel_path);

KH = kernels;
numsample = size(KH,3);

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
name_dropped = strcat(result_dir + "/dropped_kernels");
csvwrite(name_dropped, dropped)

K = 20; %number of neighbors, usually (10~30)
T = 20; %number of iterations, usually (10~20)

snf_kmeans_loc = result_dir + "/pamogk-snf-kmeans-k=" + int2str(numclass);
snf_spectral_loc = result_dir + "/pamogk-snf-spectral-k=" + int2str(numclass);

W = SNF(KH_SNF_cell, K, T);
[H_normalized_snf] = mykernelkmeans(W, numclass);
km_snf = normalized_kmeans(H_normalized_snf, numclass);
csvwrite(snf_kmeans_loc,km_snf)

group = SpectralClustering(W, numclass);
csvwrite(snf_spectral_loc,group)

kmeans_name = result_dir+"/pamogk-kmeans-k="+int2str(numclass);

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
%kmeans_name = result_dir+"/pamogk-kmeans-k="+int2str(numclass);
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
    name = strcat(result_dir + "/pamogk-mkkm-k="+int2str(numclass)+"-loglambda="+ num2str(lambdasetname(il)));
    csvwrite(name,indx2)
    name_gamma = strcat(result_dir + "/pamogk-mkkm-k="+int2str(numclass)+"-loglambda="+ num2str(lambdasetname(il))+"-weights");
    gammas = [stayed gamma2];
    csvwrite(name_gamma,gammas)
    name_obj = strcat(result_dir + "/pamogk-mkkm-k="+int2str(numclass)+"-loglambda="+ num2str(lambdasetname(il))+"-obj");
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

