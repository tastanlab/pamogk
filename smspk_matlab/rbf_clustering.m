clear
clc;
warning off

path = '';
addpath(genpath(path));
%s = RandStream('34956');
%dataName = 'ORL'; %%% flower17; proteinFold;  ucsd-mit_caltech-101-mkl
%load([path,dataName,'_Kmatrix'],'KH','Y','P');
%alpha_set = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
gamma_set = ["median2"];
for gamma = gamma_set
    folder_name = "rbf_pw_kirc_all/Experiment1-label=1-gamma="+gamma+"-norm=True";
    folderNameOut = folder_name+"/labels";
    if ~exist(folderNameOut, 'dir')
       mkdir(folderNameOut);
    end

    %KH_rs_over = readNPY(folder_name + '/rnaseq-over.npy');
    %KH_rs_under = readNPY(folder_name + '/rnaseq-under.npy');
    KH_rs_over_under = readNPY(folder_name + '/rnaseq-over-under-gamma=6.37e-03.npy');
    
    %KH_rp_over = readNPY(folder_name + '/rppa-over.npy');
    %KH_rp_under = readNPY(folder_name + '/rppa-under.npy');
    KH_rp_over_under = readNPY(folder_name + '/rppa-over-under-gamma=1.43e-01.npy');
    
    KH_som = readNPY(folder_name + '/som-gamma=4.35e-02.npy');

    numsample = size(KH_som,2);

    KH = KH_rs_over_under;
    %KH = KH_rs_over;
    %KH = cat(3,KH,KH_rs_under);
    KH = cat(3,KH,KH_rp_over_under);
    %KH = cat(3,KH,KH_rp_over);
    %KH = cat(3,KH,KH_rp_under);
    KH = cat(3,KH,KH_som);
    
    KH = kcenter(KH);
    KH = knorm(KH);
    %numclass = length(unique(Y));
    numker = size(KH,3);
    for numclass = [2,3,4,5]
        %%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gamma0 = ones(numker,1)/numker;
        avgKer = mycombFun(KH,gamma0);
        [H_normalized1] = mykernelkmeans(avgKer,numclass);
        indx = untitled(H_normalized1,numclass);
        name = folderNameOut+"/rbf-all-3ker-kmeans-"+int2str(numclass)+"lab";
        csvwrite(name,indx)

        %%%%%%%%%%---AAAI-16----%%%%%%%%
        lambdaset2 = 2.^[-15:3:15];
        %lambdaset2 = 2.^[5:2:6];

        M = calculateM(KH);
        for il =1:length(lambdaset2)
            [H_normalized2,gamma2,obj2] = myregmultikernelclustering(KH,M,numclass,lambdaset2(il));
            indx2 = untitled(H_normalized2,numclass);
            name = strcat(folderNameOut+"/rbf-all-3ker-mkkm-"+int2str(numclass)+"lab-lambda"+ num2str(il));
            csvwrite(name,indx2)
        end
    end
end

%res(:,2) = [max(accval2);max(nmival2);max(purval2)];

%% Contact 1022xinwang.liu@gmail.com
%%%%%%%-----Citation----------------------%%%%%%%
%% Xinwang Liu, Yong Dou, Jianping Yin, Lei Wang, En Zhu: 
%% Multiple Kernel k-Means Clustering with Matrix-Induced Regularization. AAAI 2016: 1888-1894

