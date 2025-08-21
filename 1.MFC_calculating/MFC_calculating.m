clc;
clear all;
NumSub = 364;
load '/1.MFC_calculating/Label_7net_5k.mat'  %7net label
load '/1.MFC_calculating/myelin.mat'  %myelin index:NumVertex*NumSub
Sublist = importdata('/data/users/wliu/demo_dHCP_Analysis/TermList_myelin.txt');  %subject information
subj = string(Sublist.textdata);
sess = string(Sublist.data);

Ind_notNuc = find(Label_7net_5k > 0);  %remove subcortical data
NumVertex = size(Ind_notNuc,1);
Ind_utri = find(triu(ones(NumVertex),1)); %matrix vectorization

%% FC
for i = 1:NumSub
    foldername = sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s',subj(i),sess(i));
    filename = sprintf('sub-%s_ses-%s_hemi-LR_BOLD_correlation.hemi_5k.dconn.nii',subj(i),sess(i)); %functional connectome
    filepath = fullfile(foldername,filename);
    sub_FC = ciftiopen(filepath);
    sub_FC = sub_FC.cdata;

    sub_FC = 0.5 * log((1+sub_FC) ./ (1-sub_FC));   %Fisher z-transform
    sub_FC=sub_FC(Ind_notNuc,Ind_notNuc);   %8589*8589
    subs_FC_Z_utri(:,i)=sub_FC_FisherZ(Ind_utri);   %matrix vectorization

    save(sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(i),sess(i)),'sub_FC');
    fprintf('%dth subject is OK~\n',i);
end

%% gMC
myelin_z_Insub = zscore(myelin,0,1);  %z-score
gMC = corr(myelin_z_Insub');  %Pearson correlation inter-subject
gMC = 0.5 * log((1+gMC) ./ (1-gMC));  %Fisher z-transform

%% sMC
for i=1:NumSub
    sub_myelin = zscore(myelin(:,i));  %z-score
    sub_MCN = sub_myelin * sub_myelin';  %Pearson correlation intra-subject
    sub_MCN = sub_MCN(Ind_utri);   %matrix vectorization
    
    save(sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(i),sess(i)),'sub_MCN');
    fprintf('%dth subject is OK~\n',i);
end

%% Vertex-level gMFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));  %Fisher-Z FC
    load(FCname);

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);   %matrix:8589*8589
        end
    end

    for m=1:length(Ind_notNuc)
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),gMC(m,[1:m-1,m+1:end])']);
        MFC_gs_v_b(m,:,k) = b;
        sub_MFC_gs_v_re(m,:) = r;
        MFC_gs_v_stats(m,:,k) = stats;
        MFC_gs_v_R2(m,k) = stats(1);
    end
    fprintf('%dth subject is OK~\n',k);
end

%% Vertex-level sMFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));
    load(FCname);
    MCNname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(k),sess(k));
    load(MCNname);

    MCN = zeros(length(Ind_notNuc),length(Ind_notNuc));
    MCN(Ind_utri) = sub_MCN;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            MCN(i,j) = MCN(j,i);  %matrix:8589*8589
        end
    end

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);   %matrix:8589*8589
        end 
    end

    for m=1:length(Ind_notNuc)
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',[ones(length(Ind_notNuc)-1,1),MCN(m,[1:m-1,m+1:end])']);
        MFC_ss_v_b(m,:,k) = b;
        sub_MFC_ss_v_re(m,:) = r;
        MFC_ss_v_stats(m,:,k) = stats;
        MFC_ss_v_R2(m,k) = stats(1);
    end
    fprintf('%dth subject is OK~\n',k);
end

%% Vertex-level MFC
for k=1:NumSub
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));
    load(FCname);
    MCNname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(k),sess(k));
    load(SCNname);

    MCN = zeros(length(Ind_notNuc),length(Ind_notNuc));
    MCN(Ind_utri) = sub_MCN;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            MCN(i,j) = MCN(j,i);  %matrix:8589*8589
        end
    end

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for i = 2:length(Ind_notNuc)
        for j = 1:i-1
            FC(i,j) = FC(j,i);  %matrix:8589*8589
        end
    end

    for m=1:length(Ind_notNuc)
        X_MCN = double([ones(length(Ind_notNuc)-1,1),zscore(gMC(m,[1:m-1,m+1:end])'),zscore(MCN(m,[1:m-1,m+1:end])')]);
        [b,~,r,~,stats]=regress(FC(m,[1:m-1,m+1:end])',X_MCN);
        MFC_gss_v_b(m,:,k) = b;
        sub_MFC_gss_v_re(m,:) = r;
        MFC_gss_v_stats(m,:,k) = stats;
        MFC_gss_v_R2(m,k) = stats(1);
    end
    fprintf('%dth subject is OK~\n',k);
end

