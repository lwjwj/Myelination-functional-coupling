clc;
clear all;
Sublist = importdata('/data/users/wliu/demo_dHCP_Analysis/TermList_myelin.txt');
subj = string(Sublist.textdata);
sess = string(Sublist.data);
load '/3.Distance_dependence_analysis/Age_364.mat'

%%%% geodesic distance of midthickness surafce 
load('/data/users/wliu/demo_dHCP_Analysis/distance.mat','distance_left')  %left sphere
load('/data/users/wliu/demo_dHCP_Analysis/distance.mat','distance_right')  %right sphere

%% distance compartment divided by number
%%%%%left sphere
for i=1:50 
    Ind_sort(i,1)=floor(2+85.8*(i-1));
    Ind_sort(i,2)=floor(2+85.8*i)-1;  %42.9
end
for i=1:size(distance_left,1)
    distance_node = distance_left(i,:);
    [~,sortedInd] = sort(distance_node,'descend');
    for j=1:50
        subnode_left{i,j} = sortedInd(Ind_sort(j,1):Ind_sort(j,2));
    end
end
%%%%%right sphere
for i=1:50
    Ind_sort(i,1)=floor(2+85.94*(i-1));
    Ind_sort(i,2)=floor(2+85.94*i)-1;    %42.97
end
for i=1:size(distance_right,1)
    distance_node = distance_right(i,:);
    [~,sortedInd] = sort(distance_node,'descend');
    for j=1:50
        subnode_right{i,j} = sortedInd(Ind_sort(j,1):Ind_sort(j,2));
    end
end

%% culculate MFC of each distance compartment
for k=1:364
    FCname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_FC.mat',subj(k),sess(k));
    load(FCname);
    MCNname=sprintf('/data/users/wliu/demo_dHCP_Analysis/sub-%s/ses-%s/sub_MCN.mat',subj(k),sess(k));
    load(MCNname);

    MCN = zeros(length(Ind_notNuc),length(Ind_notNuc));
    MCN(Ind_utri) = sub_MCN;
    for m = 2:length(Ind_notNuc)
        for n = 1:m-1
            MCN(m,n) = MCN(n,m);  %matrix
        end
    end

    FC = zeros(length(Ind_notNuc),length(Ind_notNuc));
    FC(Ind_utri) = sub_FC;
    for m = 2:length(Ind_notNuc)
        for n = 1:m-1
            FC(m,n) = FC(n,m);   %matrix
        end
    end

    for i=1:size(distance_left,1)
        for j=1:50
        subnode_left_now = subnode_left{i,j};
        MCNg_left = gMC(i,subnode_left_now)';
        MCNs_left = MCN(i,subnode_left_now)';
        FCs_left = FC(i,subnode_left_now)';
        X_MCN = double([ones(size(MCNg_left,1),1),zscore(MCNg_left),zscore(MCNs_left)]);
        MFC_gss_left(i,k,j) = regress(FCs_left,X_MCN);   %MFC of each distance compartment
        end
    end

    for i=1:size(distance_right,1)
        for j=1:50
        subnode_right_now = subnode_right{i,j};
        MCNg_right = gMC(i+4291,subnode_right_now+4291)';
        MCNs_right = MCN(i+4291,subnode_right_now+4291)';
        FCs_right = FC(i+4291,subnode_right_now+4291)';
        X_MCN = double([ones(size(MCNg_right,1),1),zscore(MCNg_right),zscore(MCNs_right)]);
        MFC_gss_right(i,k,j) = regress(FCs_right,X_MCN);    %MFC of each distance compartment
        end
    end

    fprintf('%dth sub is OK~~~~~~~~~~~~~~~~~\n',k);
end

%% average MFC across-sub
for j=1:50
MFC_gss_left_mean(:,j) = mean(MFC_gss_left(:,:,j),2);   % mean_MFC of each distance compartment
MFC_gss_right_mean(:,j) = mean(MFC_gss_right(:,:,j),2);
end

%% growth enditivity of MFC_gss_left/right_distance_7net
net = Label_7net_5k(Ind_notNuc);
for i=1:50
    for j=1:7
        for k=1:364
            MFC_gss_current_left(:,1) = MFC_gss_left(:,k,i);
            MFC_gss_current_right(:,1) = MFC_gss_right(:,k,i);
            MFC_NetMedian_indiv_50_left(i,j,k) = median(MFC_gss_current_left(find(net(1:4291) == j)),1);  %net_median MFC of each distance compartment
            MFC_NetMedian_indiv_50_right(i,j,k) = median(MFC_gss_current_right(find(net(4292:8589) == j)),1);  %net_median MFC of each distance compartment
        end
    end
end
for i=1:50
    for j=1:7
        %left sphere
        [b,~,~,~,~]=regress(MFC_NetMedian_indiv_50_left(i,j,:),[ones(364,1),PMA364]);
        b_PMA_50_left(i,j)=b(2);  %growth enditivity
        %right sphere
        [b,~,~,~,~]=regress(MFC_NetMedian_indiv_50_right(i,j,:),[ones(364,1),PMA364]);
        b_PMA_50_right(i,j)=b(2); %growth enditivity  
    end 
end

