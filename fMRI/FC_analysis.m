%% FC Analysis
%1. load in relevant data
load('colormaps.mat')
%load in schaef 400 FC
cd 'C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Graph_Theory\schaef_400\'
load('schaef_400_fc_all.mat') %schaef 400 parcellation
%create FC all matrix
cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\schaef_400'
%1000 parcellation
cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\timeseries\schaef_1000'
filename=dir('*_timeseries.mat');
%fc_all = zeros(length(filename),1102,1102); 1000 parcellation
fc_all = zeros(length(filename),502,502);

for ii=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(ii).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum)
    ses = split(2); %session ses-..
    ses = cell2mat(ses);
    load([subject_file]); %load in time-series; ts variable
    ts = ts(6:end -5,:); %remove first/last 5 time-points for noise
    %functional connectivity
    ts_corr = corr(ts); %ROI x ROI matrix
    fc_all(ii,:,:)=ts_corr; %all subjects FC into one
end



% load clinical data
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Demographic_Data/V1_MRI_all_demos.mat');
%find IDs that have 'NaN'
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
subjectV1MRIdata(sub_remove,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1);
health_sub =bin_delirium_all~=1; 
health_sub = find(health_sub==1);

fc_all(sub_remove,:,:)=[];


%% 2. Group differences
fc_del = fc_all(delirium_sub,:,:);
fc_health = fc_all(health_sub,:,:);


nROIs = size(fc_all,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_all,120,nROIs*nROIs);
fc_edges_all = flat_fc_all(:,template);
fc_edges_del = fc_edges_all(delirium_sub,:); %fc edges lower tril delirium subs
fc_edges_health = fc_edges_all(health_sub,:);
nEdges = size(fc_edges_all,2);

% run permutation across fc edges - this is the correct way! - 9/09/24
for i=1:nEdges
    [sig_fc_edges(i,:),pval_fc_edges(i,:)]=perm_code(fc_edges_del(:,i),fc_edges_health(:,i),1000);
end
% convert into matrix
sig_fc_edge_mat = matify(sig_fc_edges,502);

sig_fc_del= squeeze(mean(fc_del)).*sig_fc_edge_mat;
sig_fc_health= squeeze(mean(fc_health)).*sig_fc_edge_mat;

diff_fc_del_health = sig_fc_del - sig_fc_health;

%% 3. Variability of FC
fc_variability = std(fc_all, 0, 1); % Std deviation across subjects
fc_variability = squeeze(fc_variability); % Resulting size: 502 × 502
% variability 
fc_var_del = squeeze(std(fc_del,0,1)); % 502 x 502
fc_var_health = squeeze(std(fc_health,0,1)); 
%regional variability
region_var_del = mean(fc_var_del,2); 
region_var_health = mean(fc_var_health,2);

% figure of histogram of variability
color_del = [0.21568627450980393 0.7607843137254902 0.8509803921568627]; %blue
color_health = [0.9686274509803922 0.6039215686274509 0.4196078431372549]; %orange
figure;
set(gcf,'color','w')
hold on;
histogram(region_var_del, 'BinWidth', 0.01, 'FaceColor', color_del, 'FaceAlpha', 0.5); % Group Delirium
histogram(region_var_health, 'BinWidth', 0.01, 'FaceColor', color_health, 'FaceAlpha', 0.5); % Group Non-del
hold off;
xlabel('FC Variability per Region');
ylabel('Frequency');
legend({'Delirium', 'Non-del'});
title('Distribution of Functional Connectivity Variability Across Regions');

% subject to subject variability of FC
% Compute standard deviation across all FC values per subject
fc_subject_var = squeeze(std(fc_all, 0, [2 3])); % 120 × 1
%subject var
fc_subject_var_del = squeeze(std(fc_del, 0, [2 3]));
fc_subject_var_health = squeeze(std(fc_health, 0, [2 3]));
figure;
set(gcf,'color','w')
hold on;
histogram(fc_subject_var_del, 'BinWidth', 0.01, 'FaceColor', color_del, 'FaceAlpha', 0.5); % Group A
histogram(fc_subject_var_health, 'BinWidth', 0.01, 'FaceColor', color_health, 'FaceAlpha', 0.5); % Group B
hold off;

xlabel('FC Variability per Subject');
ylabel('Frequency');
legend({'Delirium', 'Non-Delirium'});
title('Distribution of Functional Connectivity Variability Across Subjects');

%% 4. Standard T-test with FDR Correction
[h_fc_edges2,pval_fc_edges2]=ttest2(fc_edges_del,fc_edges_health);
h_fc_edges2_mat = matify(h_fc_edges2',502);

sig_fc_del2= squeeze(mean(fc_del)).*h_fc_edges2_mat;
sig_fc_health2= squeeze(mean(fc_health)).*h_fc_edges2_mat;
diff_fc_del_health2 = sig_fc_del2 - sig_fc_health2;

%% Linear model in R - through artemis
%just get unique connections in the connectivity matrix
% only take lower triangle - unique FC edges
nROIs = size(fc_all,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_all,120,nROIs*nROIs);
fc_edges_all = flat_fc_all(:,template); %flat across fc edges only

%avg_fc_edges = mean(fc_edges_all,1);
%mat = matify(avg_fc_edges',502); %
b = 1:size(fc_edges_all,2);
c = vertcat(b,fc_edges_all); %need additional column read it into it
writematrix(c, 'fc_all_subjects.csv');
% load in the data



%% Average within-between FC across Networks
load('schaef_order.mat')
load('schaef_order8_networks.mat')

%avg. across 8 networks
order_diff_del_health = diff_fc_del_health(order_8_nets(:,2),order_8_nets(:,2));

tril_diff = tril(order_diff_del_health);
%avg within
 avg_within_diff = zeros(11,1);
 tril_diff_replace = tril_diff;
%already ordered
order_diff_del_health_replace = order_diff_del_health;
order_diff_del_health_replace(order_diff_del_health==0)=NaN;
%avg difference across and within networks
 avg_within_diff = zeros(11,11);
 tril_diff_replace = tril_diff;
 tril_diff_replace(tril_diff_replace==0)=NaN;
    for i=1:11
        for k=1:11
        avg_within_diff(i,k)= mean(tril_diff_replace(order_8_nets(:,3)==i,order_8_nets(:,3)==k),"all","omitmissing");
        end
    end




avg_within_all = zeros(120,11,11);
for k=1:120
    b = squeeze(fc_all(k,:,:));
    b_replace = b;
    b_replace(b==1)=NaN; %remove self-connections
    avg_within = zeros(11,1);
    for i=1:8
        for ii=1:8
        avg_within(i,ii)= mean(b_replace(order_8_nets(:,1)==i,order_8_nets(:,1)==ii),"all","omitmissing");
        end
     end
    avg_within_all(k,:,:)=avg_within';

end

avg_within_all_cortex = zeros(120,8,8);
for k=1:120
    b = squeeze(fc_all(k,:,:));
    b_replace = b;
    b_replace(b==1)=NaN; %remove self-connections
    avg_within = zeros(8,8);
    for i=1:8
        for ii=1:8
        avg_within(i,ii)= mean(b_replace(order_8_nets(:,1)==i,order_8_nets(:,1)==ii),"all","omitmissing");
        end
     end
    avg_within_all_cortex(k,:,:)=avg_within';

end
flat_avg_witihin_all_cortex = reshape(avg_within_all_cortex,120,8*8);
writematrix(flat_avg_within_all_cortex,'flat_avg_within_across_nets_cort.csv');


avg_within_all_del = avg_within_all(delirium_sub,:,:);
avg_within_all_health = avg_within_all(health_sub,:,:);
nNets = 11;
for i=1:nNets
   for k=1:nNets
   [sig_avg_nets(i,k),pval_avg_nets(i,k)]=perm_code(avg_within_all_del(:,i,k),avg_within_all_health(:,i,k),1000);
   end
end


%avg across network connections

% not ordered - delirium subject avg. networks
avg_net_del_sub = zeros(31,11,11);
for ss=1:size(fc_del,1) 
    a = squeeze(fc_del(ss,:,:));
    for i=1:11
        for k=1:11
        avg_within(i,k) = mean(a(order_8_nets(:,1)==i,order_8_nets(:,1)==k),'all','omitmissing'); %take avg. across entire network
        end
    end
    avg_net_del_sub(ss,:,:)=avg_within;
end
% healthy avg. networks
avg_net_health_sub = zeros(89,11,11);
for ss=1:size(fc_health,1) 
    a = squeeze(fc_health(ss,:,:));
    for i=1:11
        for k=1:11
        avg_within(i,k) = mean(a(order_8_nets(:,1)==i,order_8_nets(:,1)==k),'all','omitmissing'); %take avg. across entire network
        end
    end
    avg_net_health_sub(ss,:,:)=avg_within;
end
%significancen - across the networks
nNets = 11;
for i=1:nNets
   for k=1:nNets
   [sig_avg_nets(i,k),pval_avg_nets(i,k)]=perm_code(avg_net_del_schaef_cortex(:,i,k),avg_net_health_schaef_cortex(:,i,k),1000);
   end
end

%avg. within and across networks for delirium and healthy
sig_fc_del= squeeze(mean(fc_del)).*sig_fc_edge_mat;
sig_fc_del_order = sig_fc_del(order_8_nets(:,2),order_8_nets(:,2));
tril_sig_fc_del = tril(sig_fc_del_order);
tril_sig_fc_del(tril_sig_fc_del==0)=NaN;
nNets = 11;
for i=1:nNets
   for k=1:nNets
   avg_net_sig_del(i,k) = mean(tril_sig_fc_del(order_8_nets(:,3)==i,order_8_nets(:,3)==k),'all','omitmissing'); %take avg. across entire network
   end
end
% sig. fc health
sig_fc_health= squeeze(mean(fc_health)).*sig_fc_edge_mat;
sig_fc_health_order = sig_fc_health(order_8_nets(:,2),order_8_nets(:,2));
tril_sig_fc_health = tril(sig_fc_health_order);
tril_sig_fc_health(tril_sig_fc_health==0)=NaN;
for i=1:nNets
   for k=1:nNets
   avg_net_sig_health(i,k) = mean(tril_sig_fc_health(order_8_nets(:,3)==i,order_8_nets(:,3)==k),'all','omitmissing'); %take avg. across entire network
   end
end

%generate boxcharts for each within connection for delirium and healthy
for i=1:nNets
    within_net_del_sub(:,i)=avg_net_del_sub(:,i,i);
end
for i=1:nNets
    within_net_health_sub(:,i)=avg_net_health_sub(:,i,i);
end

figure
set(gcf,'color','w')
subplot(1,2,1)
boxchart(within_net_del_sub) %plot the distribution 
title('Avg. within net Delirium')
subplot(1,2,2)
boxchart(within_net_health_sub) %plot the distribution 
title('Avg. within net health')

A = within_net_del_sub;
B = within_net_health_sub;
% Padding the smaller dataset (A) with NaN values
max_rows = max(size(A, 1), size(B, 1));
A_padded = [A; NaN(max_rows - size(A, 1), size(A, 2))];
B_padded = [B; NaN(max_rows - size(B, 1), size(B, 2))];
% Combine datasets
combined_data = [A_padded, B_padded];  % Concatenate along columns
% Create a grouping variable for the x-axis

% Sample data: Dataset A (5 rows, 3 columns) and Dataset B (7 rows, 3 columns)
A = rand(5, 3);  % Dataset A
B = rand(7, 3);  % Dataset B

% Padding the smaller dataset (A) with NaN values
max_rows = max(size(A, 1), size(B, 1));
A_padded = [A; NaN(max_rows - size(A, 1), size(A, 2))];
B_padded = [B; NaN(max_rows - size(B, 1), size(B, 2))];

% Combine datasets
combined_data = [A_padded, B_padded];  % Concatenate along columns

% Create a grouping variable for the x-axis
grouping = [repmat(1:3, max_rows, 1), repmat(4:6, max_rows, 1)];
grouping = grouping(:);  % Flatten the grouping array

% Flatten the combined data for plotting
data = combined_data(:);

% Your data
A = rand(31, 11);  % Dataset A (31 rows, 11 columns)
B = rand(89, 11);  % Dataset B (89 rows, 11 columns)

% Padding with NaN to make the number of rows equal
max_rows = max(size(A, 1), size(B, 1));
A_padded = [A; NaN(max_rows - size(A, 1), size(A, 2))];  % 89x11 matrix
B_padded = [B; NaN(max_rows - size(B, 1), size(B, 2))];  % 89x11 matrix

% Combine the data by concatenating them column-wise
combined_data = [A_padded, B_padded];  % Now a 89x22 matrix (11 columns from A, 11 from B)

% Flatten the combined data for plotting
data = combined_data(:);

% Create a grouping variable for each column
% Create a grouping variable to assign side-by-side groups
grouping = zeros(max_rows, 22);
for col = 1:11
    % For each column in A, assign a group 2*col-1
    grouping(:, 2*col-1) = col;  % A group for this column
    % For each column in B, assign a group 2*col
    grouping(:, 2*col) = col + 0.5;  % B group for this column
end

grouping = grouping(:);  % Flatten the grouping array

% Create the boxchart
figure;
boxchart(grouping, data, 'GroupByColor', repelem([1 2], numel(grouping)/2));
xlabel('Columns');
ylabel('Values');
legend({'Dataset A', 'Dataset B'});

% above significance was not a lot - try just plotting the average fc edges
% into avg. nets




%binarize the sig. edges only - then take the average of the difference



% within network count/percentage of connections that are significantly
% either positive or negative 
order_diff_del_health
sum_pos_nets = zeros(11,1);
for i=1:nNets
    sum_pos_nets(i,:) = sum(order_diff_del_health(order_8_nets(:,3)==i,order_8_nets(:,1)==i)>0);
    %sum_neg_nets(i,:) = sum(order_diff_del_health(order_8_nets(:,3)==i,order_8_nets(:,1)==i)<0 & ~=0);

end

%data to be exported
column_labels = {'Visual','Somato-motor','Dorsal-attention','Sal-Vent attention','Limbic','Control','Default','Temp-parietal','Subcort','Cerebellum','AAS Nuclei'};


writematrix(within_net_del_sub,'avg_within_net_delirium.csv');

writematrix(within_net_health_sub,'avg_within_net_health.csv');


% avg. within permuted linear model with IL-8 interaction
filename = "permute_lm_coef_pval_interact_avg_within_cort_nets_il8.csv";
permute_coefs = readmatrix(filename);
sig_permute = double(permute_coefs(:,3)<0.05); %pvals interaction IL8
sig_coef = sig_permute.*permute_coefs(:,4);
sig_coef_reshape = reshape(sig_coef,8,8);

figure
imagesc(sig_coef_reshape)
colprmap(Color_Salmon_Blue2)
xticks([])
yticks([])
hold on
line([0,12],[1.5,1.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([1.5,1.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[2.5,2.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([2.5,2.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[3.5,3.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([3.5,3.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[4.5,4.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([4.5,4.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[5.5,5.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([5.5,5.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[7.5,7.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([7.5,7.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[9.5,9.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([9.5,9.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[10.5,10.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([10.5,10.5],[0,12],'Color','black','LineWidth',1)


%% salient-ventral attention network
load('significant_permuted_fc.mat')
load('schaef_order.mat')
load('schaef_order8_networks.mat')

diff_fc_within_salvent = diff_fc_del_health(order_8_nets(:,1)==4,order_8_nets(:,1)==4);
diff_fc_salvent = diff_fc_del_health(:,order_8_nets(:,1)==4);


[aa,bb] = sort(order_salvent(:,1));
[aa,bb] = sort(order_salvent(:,3));

%order into A/B network salvent
figure
set(gcf,'color','w')
imagesc(diff_fc_within_salvent(order_salvent(:,2),order_salvent(:,2)))
colormap(Color_Salmon_Blue2)

figure
set(gcf,'color','w')
imagesc(diff_fc_within_salvent(order_salvent(:,4),order_salvent(:,4)))
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[18.5,18.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[32.5,32.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
xticks([])
yticks([])

figure
set(gcf,'color','w')
%imagesc(diff_fc_salvent(voltron_order(1:400,3),order_salvent(:,4)))
imagesc(diff_fc_salvent(voltron_order(1:400,3),new_order_salvent(:,2))) %new order, same as chord plot
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
%line([47,47],[1,400],'Color','black','LineWidth',1)
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets

%line([0,502],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
%line([0,502],[18.5,18.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
%line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
%line([0,502],[32.5,32.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri

line([39.5,39.5],[0,502],'Color','black','LineWidth',1)
line([49.5,49.5],[0,502],'Color','black','LineWidth',1)
line([50.5,50.5],[0,502],'Color','black','LineWidth',1)

%line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
%line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
xticks([])
yticks([])

aa = mean(diff_fc_salvent(1:400,:),2);
RB_surf_schaef_pink(aa)

%need to order
order_diff_salvent = diff_fc_salvent(order_8_nets(1:400,2),new_order_salvent(:,2));
writematrix(order_diff_salvent,'order_diff_salvent_wholebrain.csv')

% create plot of DMN and Salvent Relationship
load('dmn_order_new.mat') %sub-regions of DMN
load('dmn_order_subnets.mat') %the sub-network DMN order

DMN_Salvent1 = diff_fc_salvent((order_8_nets(1:400,1)==7),new_order_salvent(:,2));

DMN_Salvent2 = DMN_Salvent1%need to reorder to DMN order

figure
set(gcf,'color','w')
imagesc(DMN_Salvent1(order_subnets,:)) %order subnetworks of DMN
colormap(Color_Salmon_Blue2)
hold on
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([39.5,39.5],[0,502],'Color','black','LineWidth',1)
line([49.5,49.5],[0,502],'Color','black','LineWidth',1)
line([50.5,50.5],[0,502],'Color','black','LineWidth',1)
%DMN subnet labels
line([0,502],[12.5,12.5],'Color','black','LineWidth',1)
line([0,502],[41.5,41.5],'Color','black','LineWidth',1)
xticks([])
yticks([])

%sub-regions of DMN
figure
set(gcf,'color','w')
imagesc(DMN_Salvent1(new_order_dmn,:)) %order subnetworks of DMN
colormap(Color_Salmon_Blue2)
hold on
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([39.5,39.5],[0,502],'Color','black','LineWidth',1)
line([49.5,49.5],[0,502],'Color','black','LineWidth',1)
line([50.5,50.5],[0,502],'Color','black','LineWidth',1)
%DMN sub-regions labels
%order DMN into regional grps
line([0,502],[9.5,9.5],'Color','black','LineWidth',1) %IPL
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %PFCd
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %PFCm
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %PFCl/v
line([0,502],[59.5,59.5],'Color','black','LineWidth',1) %pCUNPCC
line([0,502],[69.5,69.5],'Color','black','LineWidth',1) %Temp lobe
line([0,502],[74.5,74.5],'Color','black','LineWidth',1) %DefaultC_Rsp


% run permuted lm for salvent to rest of cortex
fc_all_salvent_within = fc_all(:,order_8_nets(:,1)==4,order_8_nets(:,1)==4);
fc_all_salvent_brain = fc_all(:,1:400,order_8_nets(:,1)==4);
%need to flatten fc_dmn - but remove the correlations to self
a = reshape(fc_all_salvent_within,120,51*51);
columns_to_remove = all(a == 1);
a(:, columns_to_remove) = []; %120 x10020
b = 1:size(a,2);
c = vertcat(b,a);
writematrix(c,'fc_all_salvent_wihtin.csv');
%load file
filename = "permute_lm_coef_pval_interact_salvent_within_il8.csv";
%interaction IL8 beta coefs
permute_coefs = readmatrix(filename);
sig_permute = double(permute_coefs(:,3)<0.05); %interact is less than 0.05
sig_coef = sig_permute.*permute_coefs(:,4);
% get it into the correct shape
lower_tri_vec =sig_coefs;
A= zeros(51,51);
A(tril(true(51), -1)) = lower_tri_vec;
mat_coef = A + A.';

%figure
figure
set(gcf,'color','w')
imagesc(mat_coef(order_salvent(:,4),order_salvent(:,4)))
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[18.5,18.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[32.5,32.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
xticks([])
yticks([])
%lm permuted salvent to cortex
fc_all_salvent_brain = fc_all(:,1:400,order_8_nets(:,1)==4);
a = reshape(fc_all_salvent_brain,120,400*51);
columns_to_remove = all(a == 1);
a(:, columns_to_remove) = []; %120 x10020
b = 1:size(a,2);
c = vertcat(b,a);
writematrix(c,'fc_all_salvent_brain.csv');

% Export the data for analysis in LM
tril_within_salvent = zeros(120,51,51);
for i = 1:120
    % Extract the lower triangle of the current 51x51 slice
    tril_within_salvent(i,:,:) = tril(squeeze(fc_all_salvent_within(i, :, :)), -1);  % Extracts lower triangle without the diagonal
end

flat_tril_within_salvent =reshape(tril_within_salvent,120,51*51);
% Find columns with any zeros across all slices in the 3D array
cols_with_zero = any(flat_tril_within_salvent == 0);

% Remove columns with zeros from all slices
flat_tril_within_salvent(:, cols_with_zero) = [];

writematrix(flat_tril_within_salvent,'upper_tri_salvent_brain.csv');

% Load in coeff weights - to do the plots
filename = "permute_lm_coef_pval_fc_within_salvent_dsst.csv";
permute_coefs = readmatrix(filename);
sig_permute = double(permute_coefs(:,2)<0.05);
sig_coef = sig_permute.*permute_coefs(:,1); %multiply the beta coefs with the p-value

% get it into the correct shape
mat_coef=matify(sig_coef,51);

%sig. beta coefs from permute lm fc salvent within fc to tmta
figure
set(gcf,'color','w')
imagesc(mat_coef(order_salvent(:,4),order_salvent(:,4)))
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[18.5,18.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([18.5,18.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[25.5,25.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([25.5,25.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[32.5,32.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([32.5,32.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[37.5,37.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([37.5,37.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[47.5,47.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([47.5,47.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
xticks([])
yticks([])


%% Dorsal Attention Network:
fc_all_dorsatten_within = fc_all(:,[60:85,259:284],[60:85,259:284]);
% Export the data for analysis in LM
tril_within_dorsatten = zeros(120,52,52);
for i = 1:120
    % Extract the lower triangle of the current 51x51 slice
    tril_within_dorsatten(i,:,:) = tril(squeeze(fc_all_dorsatten_within(i, :, :)), -1);  % Extracts lower triangle without the diagonal
end

flat_tril_within_dorsatten =reshape(tril_within_dorsatten,120,52*52);
% Find columns with any zeros across all slices in the 3D array
cols_with_zero_dorsatten = any(flat_tril_within_dorsatten == 0);

% Remove columns with zeros from all slices
flat_tril_within_dorsatten(:, cols_with_zero_dorsatten) = [];

writematrix(flat_tril_within_dorsatten,'upper_tri_dorsatten_brain.csv');

% Load in coeff weights - to do the plots
filename = "permute_lm_coef_pval_fc_within_dorsatten_dsst.csv";
permute_coefs = readmatrix(filename);
sig_permute = double(permute_coefs(:,2)<0.05);
sig_coef = sig_permute.*permute_coefs(:,1); %multiply the beta coefs with the p-value

% get it into the correct shape
mat_coef=matify(sig_coef,52);

%sig. beta coefs from permute lm fc dorsatten within fc to tmta
figure
set(gcf,'color','w')
imagesc(mat_coef(order_dorsatten(:,3),order_dorsatten(:,3)))
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[13.5,13.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([13.5,13.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[28.5,28.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([28.5,28.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[45.5,45.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([45.5,45.5],[0,502],'Color','black','LineWidth',1) %vis Cent + Peri
xticks([])
yticks([])


%% Visual Network Analysis

%whole brain to visual network ordered
diff_vis = diff_fc_del_health(:,order_8_nets(:,1)==1);

%4 grouped subnetworks of visual network - column 2 is reorder, column 3 is
%order
diff_order_vis = diff_vis(:,order_vis(:,2));

figure
set(gcf,'color','w')
imagesc(diff_order_vis(voltron_order(1:400,3),:)) %new order, same as chord plot
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
%line([47,47],[1,400],'Color','black','LineWidth',1)
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
%vis lines
line([24.5,24.5],[0,502],'Color','black','LineWidth',1) %VisCentExtStr
line([33.5,33.5],[0,502],'Color','black','LineWidth',1) %VisPeri_ExtrSup
line([43.5,43.5],[0,502],'Color','black','LineWidth',1) %VisPeri_ExtrInf
xticks([])
yticks([])

%takes the significant avg only - removes non-significant
bb = diff_order_vis;
bb(bb==0)=NaN;
avg_aa = mean(bb(1:400,:),2,'omitmissing');
avg_aa(isnan(avg_aa))=0;

aa = mean(diff_order_vis(1:400,:),2);
RB_surf_schaef_pink(aa)

%need to order
order_diff_vis = diff_vis(order_8_nets(1:400,2),:);
writematrix(order_diff_vis,'order_diff_visual_wholebrain.csv')

%lm permuted visual to cortex
fc_all_vis_brain = fc_all(:,1:400,order_8_nets(:,1)==1);
a = reshape(fc_all_vis_brain,120,400*47);
columns_to_remove = all(a == 1);
a(:, columns_to_remove) = []; %120 x10020
b = 1:size(a,2);
c = vertcat(b,a);
writematrix(c,'fc_all_vis_brain.csv');

%load in variable
%load file
filename = "permute_lm_coef_pval_interact_vis_brain_il8.csv";
%interaction IL8 beta coefs
permute_coefs = readmatrix(filename);
sig_permute = double(permute_coefs(:,3)<0.05); %interact is less than 0.05
sig_coef = sig_permute.*permute_coefs(:,4);
%add the region back self connections
locs_self = find(columns_to_remove==1);
    numNewElements = length(locs_self);
    B = zeros(length(sig_coef)+numNewElements,1);
    originalIndex =1;
    newIndex = 1;
    for i = 1:length(sig_coef) + numNewElements
        if ismember(i, locs_self)
            B(i) = 0;  % Insert zero at the specified position
        else
            B(i) = sig_coef(originalIndex);  % Copy the original element
            originalIndex = originalIndex + 1;
        end
    end
     %sig. beta coefs from permute lm withn salvent fc to brain
  mat_coef = reshape(B,400,47);
%fig beta coef interact IL8
figure
set(gcf,'color','w')
imagesc(mat_coef(voltron_order(1:400,3),order_vis(:,2))) %new order, same as chord plot
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
%line([47,47],[1,400],'Color','black','LineWidth',1)
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
%vis lines
line([24.5,24.5],[0,502],'Color','black','LineWidth',1) %VisCentExtStr
line([33.5,33.5],[0,502],'Color','black','LineWidth',1) %VisPeri_ExtrSup
line([43.5,43.5],[0,502],'Color','black','LineWidth',1) %VisPeri_ExtrInf
xticks([])
yticks([])

%% Data for Chord Plots

%difference in connections for within plotted into a chord plot

%order - into 11 network
writematrix(order_diff_del_health,'order_sig_diff_del_health.csv');
writematrix(order_8_nets,'grouping_11_nets_fc.csv');

%Data for Salvent to rest of cortex - chord plot:



%% Reduce Dimensionality -
fc_400 = fc_all(:,1:400,1:400);
nROIs = size(fc_400,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_400,120,nROIs*nROIs);
fc_edges_all_400 = flat_fc_all(:,template); %flat across fc edges only

writematrix(fc_edges_all_400, 'flat_fc_edges_400_schaef_all_subjects.csv');

%avg within networks for cortex - flatten
flat_avg_witihin_cortex = reshape(avg_within_all_cortex,120,8*8);
writematrix(flat_avg_witihin_cortex,'avg_within_across_cort_networks.csv');
%% Create Final Figures

diff_fc_del_health = sig_fc_del - sig_fc_health;

figure
set(gcf,'color','w')
%subplot(1,2,1)
imagesc(diff_fc_del_health(voltron_order(1:400,3),voltron_order(1:400,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Avg. diff del - health')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,400],'Color','black','LineWidth',1)
line([1,400],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([1,400],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([1,400],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([1,400],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([1,400],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([1,400],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue2)

%subcortical
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
%rest to 502 are asending arousal
colormap(Color_Salmon_Blue)

%just plot the lower triangle 
order_diff_del_heallth = diff_fc_del_health(voltron_order(:,3),voltron_order(:,3));

tril_diff = tril(order_diff_del_heallth);

figure
set(gcf,'Color','w')
imagesc(tril_diff(1:400,1:400))
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,400],'Color','black','LineWidth',1)
line([1,400],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,400],'Color','black','LineWidth',1) %somat mot a/b
line([1,400],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,400],'Color','black','LineWidth',1) %dorsatten
line([1,400],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,400],'Color','black','LineWidth',1) %sal ventatten
line([1,400],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,400],'Color','black','LineWidth',1) % limbic
line([1,400],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,400],'Color','black','LineWidth',1) %all Cont A-C
line([1,400],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,400],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue)
xticklabels([])
yticklabels([])

% create an avg. plot of the connection for each network

figure
imagesc(avg_within_diff(1:8,1:8))
xticks([])
yticks([])
title('Avg. Within & Across Network')
hold on
line([0,12],[1.5,1.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([1.5,1.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[2.5,2.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([2.5,2.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[3.5,3.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([3.5,3.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[4.5,4.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([4.5,4.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[5.5,5.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([5.5,5.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[7.5,7.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([7.5,7.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[9.5,9.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([9.5,9.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[10.5,10.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([10.5,10.5],[0,12],'Color','black','LineWidth',1)
colormap(Color_Salmon_Blue2)
colorbar
%plot of grp avg. within and across networks
figure
set(gcf,'color','w')
subplot(1,3,1)
avg_net_sig_del(isnan(avg_net_sig_del))=0;
imagesc(avg_net_sig_del(1:8,1:8))
xticks([])
yticks([])
title('Avg. Within & Across Network Delirium')
hold on
line([0,12],[1.5,1.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([1.5,1.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[2.5,2.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([2.5,2.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[3.5,3.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([3.5,3.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[4.5,4.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([4.5,4.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[5.5,5.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([5.5,5.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[7.5,7.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([7.5,7.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[9.5,9.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([9.5,9.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[10.5,10.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([10.5,10.5],[0,12],'Color','black','LineWidth',1)
colormap(Color_Salmon_Blue2)
colorbar
subplot(1,3,2)
avg_net_sig_health(isnan(avg_net_sig_health))=0;
imagesc(avg_net_sig_health(1:8,1:8))
xticks([])
yticks([])
title('Avg. Within & Across Network Health')
hold on
line([0,12],[1.5,1.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([1.5,1.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[2.5,2.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([2.5,2.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[3.5,3.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([3.5,3.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[4.5,4.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([4.5,4.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[5.5,5.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([5.5,5.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[7.5,7.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([7.5,7.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[9.5,9.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([9.5,9.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[10.5,10.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([10.5,10.5],[0,12],'Color','black','LineWidth',1)
colormap(Color_Salmon_Blue2)
colorbar
subplot(1,3,3)
imagesc(avg_within_diff(1:8,1:8))
xticks([])
yticks([])
title('Avg. Difference Within & Across Network Del - Health')
hold on
line([0,12],[1.5,1.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([1.5,1.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[2.5,2.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([2.5,2.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[3.5,3.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([3.5,3.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[4.5,4.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([4.5,4.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[5.5,5.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([5.5,5.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[6.5,6.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([6.5,6.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[7.5,7.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([7.5,7.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[8.5,8.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([8.5,8.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[9.5,9.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([9.5,9.5],[0,12],'Color','black','LineWidth',1)
line([0,12],[10.5,10.5],'Color','black','LineWidth',1) %vis Cent + Peri
line([10.5,10.5],[0,12],'Color','black','LineWidth',1)
colormap(Color_Salmon_Blue2)
colorbar



%% Figures
%load in network connectivity order
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\PhD\Code\Parcellations\voltron_ordered.mat')

figure
set(gcf,'Color','w')
subplot(1,2,1)
mean_fc_del = squeeze(mean(fc_del,1));
tril_del = tril(mean_fc_del(voltron_order(:,3),voltron_order(:,3)));
imagesc(tril_del(1:400,1:400))
title('Avg. FC Delirium')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue2)
xticks([])
yticks([])
colorbar
subplot(1,2,2)
mean_fc_health = squeeze(mean(fc_health,1));
tril_health = tril(mean_fc_health(voltron_order(:,3),voltron_order(:,3)));
imagesc(tril_health(1:400,1:400))
title('Avg. FC Healthy')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
colormap(Color_Salmon_Blue2)
xticks([])
yticks([])
colorbar





figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(sig_fc_del(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Avg. FC delirium')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
%rest to 502 are asending arousal
subplot(1,2,2)
imagesc(sig_fc_health(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Avg. FC healthy')
colormap(CustomColormap)
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum


figure
set(gcf,'color','w')
%subplot(1,2,1)
imagesc(diff_fc_health_del(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Avg. diff health - del')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
%rest to 502 are asending arousal
colormap(Color_Salmon_Blue)
%plot correlation with post-op
figure
set(gcf,'color','w')
%imagesc(reshape_rho_fc_peak_drs(voltron_order(:,3),voltron_order(:,3)))
imagesc(reshape_sig_permute_corr(voltron_order(:,3),voltron_order(:,3)))
xlabel('ROIs')
ylabel('ROIs')
title('Sig. Avg. FC per edge corr with postop drs peak')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(Color_Salmon_Blue)

figure
set(gcf,'color','w')
imagesc(sig_reshape_rho_fc_peak_drs(voltron_order(:,3),voltron_order(:,3)))
xlabel('ROIs')
ylabel('ROIs')
title('Sig. Avg. FC per edge corr with postop drs peak')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(CustomColormap)


% GLM results - only sig. betas
figure
set(gcf,'color','w')
imagesc(mat_sig_beta_coeff_ROI_edges(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('Sig. Beta Coeff FC')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(CustomColormap)

% just plot negative beta coeffs -


% check within DMN plot -

figure
imagesc(diff_fc_del_health(order_8_nets(:,1)==7,order_8_nets(:,1)==7))
colormap(Color_Salmon_Blue2)