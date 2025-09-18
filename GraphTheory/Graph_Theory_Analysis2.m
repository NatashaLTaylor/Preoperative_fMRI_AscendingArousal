%% Graph Theory Analysis 2

%load in relevant data
%load in clinical data
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Demographic_Data\V1_MRI_all_demos.mat');
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
sub_remove(3,:)=[]; % need to remove 59, as that has already been removed
sub_remove_clinical = [42,43,59,63,74,93,103,104];
subjectV1MRIdata(sub_remove_clinical,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1);
health_sub =bin_delirium_all~=1; 
health_sub = find(health_sub==1);

%% dFC - permute difference


load("mtd_all_flat.mat")





%% Cartographic Profile
load('C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\Graph_Theory\schaef_400\cartographic_all.mat')
% size 125 x 101 x 101 x 210
cartograph = cartographic_all;
cartograph([sub_remove],:,:,:)=[]; %size 119 x 101 x 101 x 210
%group cartographs
del_cartograph = cartograph(delirium_sub,:,:,:);
health_cartograph = cartograph(health_sub,:,:,:);

mean_del_cart = squeeze(mean(del_cartograph,4));
mean_health_cart = squeeze(mean(health_cartograph,4));



figure
set(gcf,'Color','w')
subplot(1,2,1)
imagesc(squeeze(mean(mean_del_cart,1)))
title('Delirium Cartograph Avg. Time')
xlabel('Participation Coeff')
ylabel('Module Degree')
subplot(1,2,2)
imagesc(squeeze(mean(mean_health_cart,1)))
title('Health Cartograph Avg. Time')
xlabel('Participation Coeff')
ylabel('Module Degree')
% Diff Fig
load('colormaps2.mat')

figure
imagesc(squeeze(mean(mean_del_cart,1)-mean(mean_health_cart,1)))
 ylabel('Module Degree')
xlabel('Participation Coeff')
set(gcf,'Color','w')
title('Mean Diff Cartograph Del - Health')
colormap('Colormap_Salmon_Blue')


diff2 = mean(mean_del_cart,1)-mean(mean_health_cart,1);

% percentage of time







%% K-means clustering of the states - 
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Graph_Theory\schaef_400\participation_coeff_all.mat')
part_all(sub_remove,:,:)=[]; %remove subjects without clinical data
part_del = part_all(delirium_sub,:,:);
part_healthy = part_all(health_sub,:,:);

%reorder to export for sklearn in python
% attempt run cluster with sub x time x ROI
reorder_part_all = permute(part_all,[1,3,2]);
writematrix(reorder_part_all,"reorder_pc_all.csv")

%kmeans clustering on the pc values - for each ROI
a = squeeze(part_all(1,:,:));

idx = zeros(210,502);
C = zeros(2,502);
for i=1:502
    b= a(i,:);
    [idx(:,i),C(2,i)]=kmeans(b',2);
end


[idx,C,sumd]=kmeans(a',2);

silhouette(a',idx)

% test on 1 subject
a = squeeze(part_all(1,:,:));
% testing the number of clusters - elbow point approach
wcss = zeros(100,1);
for k = 1:100
    [~, ~, sumd] = kmeans(a', k);
    wcss(k) = sum(sumd);
end
% Plot the WCSS to find the elbow
figure;
plot(1:100, wcss, '-o');
xlabel('Number of clusters (k)');
ylabel('Within-cluster sum of squares (WCSS)');
title('Elbow Method for Optimal k');

% silhouette approach for optimising
for k = 1:10
    idx = kmeans(a', k);
    s = silhouette(a', idx);
    silhouetteAvg(k) = mean(s);
end

% Plot the silhouette values to find the optimal number of clusters
%pick the number of clusters that has the highest avg. silhouette value
figure;
plot(2:10, silhouetteAvg(2:end), '-o');
xlabel('Number of clusters (k)');
ylabel('Average Silhouette Value');
title('Silhouette Method for Optimal k');

%attempt clustering all subjects - then determine overall best number of
%clusters
silhouetteAvg_sub = zeros(119,10);
for i=1:119
    a = squeeze(part_all(i,:,:)); %load in sub pc
        % silhouette approach for optimising
    for k = 1:10
        idx = kmeans(a', k);
        s = silhouette(a', idx);
        silhouetteAvg(k) = mean(s);
    end
        silhouetteAvg_sub(i,:)=silhouetteAvg;
end
%figure plot
figure;
plot(2:10, squeeze(silhouetteAvg_sub(i,2:end)), '-o');
xlabel('Number of clusters (k)');
ylabel('Average Silhouette Value');
title('Silhouette Method for Optimal k');

figure
plot(1:119,squeeze(silhouetteAvg_sub(:,2)),'-o')
xlabel('Subjects')
ylabel('Avg. Silhouette for K=2')

% Run Clustering on Groups -
%run cluster
idx_k = zeros(210,20);
C_k = zeros()
for k = 1:20
    [idx, C] = kmeans(concat_part_del, k);
    idx_k(:,k)=idx;
    C_k()
    s = silhouette(concat_part_del, idx);
    silhouetteAvg(k) = mean(s);
end

figure;
plot(2:20, silhouetteAvg(2:end), '-o');
xlabel('Number of clusters (k)');
ylabel('Average Silhouette Value');
title('Silhouette Method for Optimal k Delirium Grp');

% run cluster on del grp
aa = permute(part_del,[3,2,1]);
reshapedData = reshape(aa,[],502*31);

% attempt run cluster with sub x time x ROI
aa = permute(part_del,[1,3,2]);
reshapedData = reshape(aa,[],502); %6510 x 502


[idx, C] = kmeans(reshapedData,3);

stateassignment = reshape(idx,31,210);
%plot for each state
for subj = 1:size(part_del,1)
    figure;
    imagesc(stateassignment(subj, :));
    colormap('jet');
    colorbar;

    % Customize the plot
    title(['State Assignment for Subject ', num2str(subj)]);
    xlabel('Time Points');
    ylabel('State');
end
% Plot the cluster centroids
figure;
plot(C');
legend(arrayfun(@(x) ['State ' num2str(x)], 1:k, 'UniformOutput', false));
title('Cluster Centroids');
xlabel('Region');
ylabel('Centroid Value');






%% Running Bin Grouped Stats

% group-wise permutation significant of PC edges
nROIs = size(part_all,2);
for i=1:nROIs
   for k=1:size(part_all,3)
   [sig_pc_time(i,k),pval_pc_time(i,k)]=perm_code(part_del(:,i,k),part_healthy(:,i,k),1000);
   end
end


mean_part_del= squeeze(mean(part_del,1));
mean_part_health= squeeze(mean(part_healthy,1));

load('schaef_order.mat')
load('colormaps2.mat')
figure
set(gcf,'color','w')
subplot(1,3,1)
a= mean_part_del.*sig_pc_time;
imagesc(a(voltron_order(:,3),:))
title('Avg. Delirium PC')
hold on
line([0,210],[47,47],'Color','black','LineWidth',1)
line([0,210],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,210],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,210],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,210],[244,244],'Color','black','LineWidth',1) %limbic
line([0,210],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,210],[384,384],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
yticklabels([])
xlabel('Time')
subplot(1,3,2)
b= mean_part_health.*sig_pc_time;
imagesc(b(voltron_order(:,3),:))
title('Avg. Health PC')
hold on
line([0,210],[47,47],'Color','black','LineWidth',1)
line([0,210],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,210],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,210],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,210],[244,244],'Color','black','LineWidth',1) %limbic
line([0,210],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,210],[384,384],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
yticklabels([])
xlabel('Time')
subplot(1,3,3)
c= a - b;
imagesc(c(voltron_order(:,3),:))
title('Avg. Diff PC Del - Health')
hold on
line([0,210],[47,47],'Color','black','LineWidth',1)
line([0,210],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,210],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,210],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,210],[244,244],'Color','black','LineWidth',1) %limbic
line([0,210],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,210],[384,384],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
yticklabels([])
xlabel('Time')


%sig different regions over time for pc and the diff



%plot onto brain
mean_diff_del_health_pc = mean(c,2);
RB_surf_schaef_pink(mean_diff_del_health_pc(1:400,1),'~') %plot the diff in PC 
RB_surf_schaef_color(mean_diff_del_health_pc(1:400,1),'~',[196 64 48],[12 125 121])

writematrix(mean_diff_del_health_pc,'mean_diff_del_health_pc_avg_time.csv')
% subject level - significant PC based on time for each ROI
for i=1:length(delirium_sub)
    sig_pc_time_del(i,:,:) = sig_pc_time.*squeeze(part_del(i,:,:));
end
sig_avg_pc_del = mean(sig_pc_time_del,3);
sig_avg_pc_time_del = mean(sig_pc_time_del,1);
for i=1:length(health_sub)
    sig_pc_time_health(i,:,:) = sig_pc_time.*squeeze(part_healthy(i,:,:));
end
sig_avg_pc_health = mean(sig_pc_time_health,3);
sig_avg_pc_time_health = mean(sig_pc_time_health,1);
writematrix(sig_avg_pc_del,'sig_avg_pc_del.csv'); %avg. across time sub x pc
writematrix(sig_avg_pc_health,'sig_avg_pc_nondel.csv');

diff_avg_pc_del_health = sig_avg_pc_time_del - sig_avg_pc_time_health;

mean_diff_avg_pc_del_health = mean(diff_avg_pc_del_health,2);
RB_surf_schaef_color(mean_diff_avg_pc_del_health(1:400,1),'~',[196 64 48],[12 125 121])








%create avg. of PC across time - 
avg_pc_del = mean(part_del,3);
avg_pc_health = mean(part_healthy,3);

% permuted difference in just the avg. PC
% run permutation across fc edges - this is the correct way! - 9/09/24
nEdges = 502;
for i=1:nEdges
    [sig_avg_pc(i,:),pval_avg_pc(i,:)]=perm_code(avg_pc_del(:,i),avg_pc_health(:,i),1000);
end

mean_diff_avg_pc = mean(avg_pc_del,1) - mean(avg_pc_health,1);

diff_avg_pc_sig = 


%output .csv for the significant avg. pc ROIs and the avg for the del and
%health pops
writematrix(sig_avg_pc,'significant_rois_avg_pc.csv');
writematrix(avg_pc_del,'avg_pc_del.csv'); %avg. across time sub x pc
writematrix(avg_pc_health,'avg_pc_nondel.csv');
%all subs together - avg pc overtime
avg_pc_all = mean(part_all,3);
writematrix(avg_pc_all,'avg_pc_all_subs.csv');
%saving output of pc x time significance
writematrix(part_del,'pc_delirium.csv'); %sub x roi x time
writematrix(part_healthy,'pc_nondel.csv'); %sub x roi x time
writematrix(sig_pc_time,'sig_permute_pc_time.csv'); %sig permute pc x time

writematrix(sig_avg_pc,'sig_avg_pc.csv'); %avg. pc over time


%% Modularity
load('modularity_all.mat')
modularity_all(sub_remove,:)=[];

modularity_del=modularity_all(delirium_sub,:);
modularity_health=modularity_all(health_sub,:);

mean_mod_del = mean(modularity_del);
mean_mod_health = mean(modularity_health);

%relationship between modularity and delirium severity

mean_modularity_all = mean(modularity_all,2);
%relationship between the two
peak_drs_postop =table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"overall_peak_drs_calc2"));
[rho_modularity,pval_modularity]=corr(mean_modularity_all,peak_drs_postop,'rows','complete');



%% Flexibility
load('flexibility_all.mat')

flex_all(sub_remove,:)=[];

[rho_flex,pval_flex]=corr(flex_all,peak_drs_postop,'rows','complete');
%29 regions are correlated with flexibility

sig_flex = double(pval_flex<0.05);
sig_rho_flex = rho_flex.*sig_flex;

locs = find(sig_flex==1);



RB_surf_schaef_pink(sig_rho_flex(1:400,1),'~')




%% Network-level Graph Theory metrics

