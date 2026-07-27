%% REDO ANALYSIS - NAT COMMS

% remove relevant subjects - MRI scanner acquisition:
% load in the following
%readmatrix('mri_scanner_current.csv');
load('/Users/ntaylor/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Demographic_Data/V1_MRI_all_demos.mat');

%sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
%sub_remove(3,:)=[]; % need to remove 59, as that has already been removed
sub_remove_clinical = [42,43,59,63,74,93,103,104];
subjectV1MRIdata_graphtheory_orig=subjectV1MRIdata;
subjectV1MRIdata_graphtheory_orig(sub_remove_clinical,:)=[];
% original indices for clinical groups
bin_delirium_all=table2array(subjectV1MRIdata_graphtheory_orig(1:size(subjectV1MRIdata_graphtheory_orig,1),"bin_delirium"));
delirium_sub_orig = find(bin_delirium_all==1);
health_sub_orig =bin_delirium_all~=1; 
health_sub_orig = find(health_sub_orig==1);

% remove 1.5T data
a= find(subjectV1MRIdata_graphtheory_orig.sub_id==49);
b= find(subjectV1MRIdata_graphtheory_orig.sub_id==50);
sub_remove_grp1=[a,b];
subjectV1MRIdata_graphtheory_update = subjectV1MRIdata_graphtheory_orig;
subjectV1MRIdata_graphtheory_update(sub_remove_grp1,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata_graphtheory_update(1:size(subjectV1MRIdata_graphtheory_update,1),"bin_delirium"));
delirium_sub_update = find(bin_delirium_all==1);
health_sub_update =bin_delirium_all~=1; 
health_sub_update = find(health_sub_update==1);

writematrix(bin_delirium_all,'delirium_idx.csv')

subject_id_updated = subjectV1MRIdata_graphtheory_update.sub_id;
% Convert numbers to zero-padded strings with "sub-" prefix
subject_id_updated = compose("sub-%03d", subject_id_updated);

writematrix(subject_id_updated,'subject_ids.txt')

% look at whether there's a significant difference for LC & NBM
% connectivity between groups
% scanner 1.5T (group 1) - sub-049 & sub-050




%----------------------------------------%
%% GRAPH THEORY REDO
%----------------------------------------%

%----------------------------------------%
%% RECALCULATE PC
%----------------------------------------%
output_dir = "/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/csf_regress_redo";
% subject-wise PC values
%% load data MTD
cd(output_dir)
filename=dir('*mtd_flat.mat');
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    split = strsplit(subject_file,'csf4thvent');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum);
    %ses = split(2); %session ses-..
    %ses = cell2mat(ses);
    
    load([subject_file]); %load in mtd
    %need to matify the flatten dFC

    for t = 1:n_timepoints
        A = matify(mtd_flat(:, t), n_ROIs);  % Reconstruct full 502x502 matrix
        mtd(:, :, t) = A;  % Only keep rows from selected regions
    end

    
%data size is 502 x 502 x 210
    gamma = 1;
    beta = 0.75;

    [ci,q,part,~,hc,~] = integration_plus5(mtd,gamma,beta);

    save([subnum '_csf4thvent_part.mat'],"part");
    %save([subnum '_' ses '_flex.mat'],"f");
    %save([subnum '_' ses '_modz.mat'],"modz");
    %save([subnum '_' ses '_modularity.mat'],"q");
    %save([subnum '_' ses '_community.mat'],"ci");
    save([subnum '_csf4thvent_cartographic.mat'],"hc");
    fprintf('Finished subject %d of %d\n', i, n_subjects);
end


% remove the 1.5T subjects
% load in all pc
% attempt run cluster with sub x time x ROI
%reorder_part_all = permute(part_all,[1,3,2]);
%writematrix(reorder_part_all,"reorder_pc_all.csv")

filename=dir('*part.mat');
%remove shorter subject windows - sub-107 & sub-164 - already not in the
%analysis
part_all = zeros(length(filename),502,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

part_all(i,:,:)= part;
end
% remove two subjects - 1.5T subjects
part_all(sub_remove_grp1,:,:)=[];
% needs to be sub x time x ROI
reorder_part_all = permute(part_all,[1,3,2]);
writematrix(reorder_part_all,'reorder_pc_all.csv')


part_del = part_all(delirium_sub_update,:,:);
part_healthy = part_all(health_sub_update,:,:);

% permutation testing for 
% group-wise permutation significant of PC edges

nROIs = size(part_all,2);
for i=1:nROIs
   for k=1:size(part_all,3)
   [sig_pc_time(i,k),pval_pc_time(i,k)]=perm_code(part_del(:,i,k),part_healthy(:,i,k),1000);
   end
end

% original files: mean_diff_del_health_pc_avg_time.csv




% subject level - significant PC based on time for each ROI
for i=1:length(delirium_sub_update)
    sig_pc_time_del(i,:,:) = sig_pc_time.*squeeze(part_del(i,:,:));
end
sig_avg_pc_del = mean(sig_pc_time_del,3);
sig_avg_pc_time_del = mean(sig_pc_time_del,1);
for i=1:length(health_sub_update)
    sig_pc_time_health(i,:,:) = sig_pc_time.*squeeze(part_healthy(i,:,:));
end
sig_avg_pc_health = mean(sig_pc_time_health,3);
sig_avg_pc_time_health = mean(sig_pc_time_health,1);
writematrix(sig_avg_pc_del,'sig_avg_pc_del.csv'); %avg. across time sub x pc
writematrix(sig_avg_pc_health,'sig_avg_pc_nondel.csv');

diff_avg_pc_del_health = sig_avg_pc_time_del - sig_avg_pc_time_health;

mean_diff_avg_pc_del_health = mean(diff_avg_pc_del_health,2);
writematrix(mean_diff_avg_pc_del_healt,'mean_diff_del_health_pc_avg_time.csv')

%all subs together - avg pc overtime
avg_pc_all = mean(part_all,3);
writematrix(avg_pc_all,'avg_pc_all_subs.csv');

% average in PC between groups - for the significantly different timepoints

%look at the sig - pc timepoint for each and take an average
mean_pc_oversigtime_del = mean(sig_avg_pc_time_del(sig_avg_pc_time_del ~=0));
mean_pc_oversigtime_health = mean(sig_avg_pc_time_health(sig_avg_pc_time_health ~=0));

% avg pc significant difference
% run permutation across fc edges - this is the correct way! - 9/09/24
avg_pc_del = mean(part_del,3);
avg_pc_health = mean(part_healthy,3);

nEdges = 502;
for i=1:nEdges
    [sig_avg_pc(i,:),pval_avg_pc(i,:)]=perm_code(avg_pc_del(:,i),avg_pc_health(:,i),1000);
end

mean_diff_avg_pc = mean(avg_pc_del,1) - mean(avg_pc_health,1);

%output .csv for the significant avg. pc ROIs and the avg for the del and
%health pops
writematrix(sig_avg_pc,'significant_rois_avg_pc.csv');
writematrix(avg_pc_del,'avg_pc_del.csv'); %avg. across time sub x pc
writematrix(avg_pc_health,'avg_pc_nondel.csv');


% reshape accordingly 
[n_subs, n_regions, n_time] = size(part_all);

% Stack: all subjects at time 1, then all subjects at time 2, etc.
pc_stacked = [];
for t = 1:n_time
    pc_stacked = [pc_stacked; part_all(:, :, t)];
end

% pc_stacked is now (n_subs * n_time) x n_regions
csvwrite('pc_data_time.csv', pc_stacked);


% effect size: 
% Compute along the subject dimension (dim = 1)
n1 = size(part_del, 1);
n2 = size(part_healthy, 1);
n_rois = size(part_del, 2);
n_time = size(part_healthy, 3);

mean1 = squeeze(mean(part_del, 1));  % [n_rois x n_time]
mean2 = squeeze(mean(part_healthy, 1));  % [n_rois x n_time]

var1 = squeeze(var(part_del, 0, 1)); % unbiased variance
var2 = squeeze(var(part_healthy, 0, 1));

pooled_sd = sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1 + n2 - 2));

cohens_d = (mean1 - mean2) ./ pooled_sd;

% Apply Hedge's g correction
correction = 1 - (3 / (4*(n1 + n2) - 9));
hedges_g = cohens_d * correction;
% visualize
sig_effect_sizes = hedges_g .* sig_pc_time;

% Visualize
imagesc(sig_effect_sizes);
colorbar;
xlabel('Time');
ylabel('ROI');
title('Hedge''s g (significant clusters only)');


% Average participation coefficient across time for each subject x ROI
group1_avg = squeeze(mean(part_del, 3));  % [n_subjects x n_rois]
group2_avg = squeeze(mean(part_healthy, 3));  % [n_subjects x n_rois]

% Then compute Cohen's d per ROI
for roi = 1:n_rois
    g1 = group1_avg(:, roi);
    g2 = group2_avg(:, roi);
    pooled_sd = sqrt(((n1-1)*var(g1) + (n2-1)*var(g2)) / (n1 + n2 - 2));
    d_per_roi(roi) = (mean(g1) - mean(g2)) / pooled_sd;
end

%----------------------------------------%
%% RECALCULATE dFC - LC & nbM
%----------------------------------------%
output_dir = "/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/csf_regress_redo";
cd(output_dir)
load('mtd_flat_all_subjects_csf4thvent.mat')
mtd_flat_csf4thvent = mtd_flat_all;
% remove 1.5T subjects
mtd_flat_csf4thvent_remove = mtd_flat_csf4thvent;
mtd_flat_csf4thvent(sub_remove_grp1,:,:)=[];

% compare the difference - old dFC & new dFC
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/AAS_MTD_permuted.mat');

data = mtd_flat_csf4thvent;  % sub x ROI x ROI x time
[n_subjects, ~, n_timepoints] = size(data);
n_ROIs = 502;
%selected_regions = [485, 486, 489, 490, 491, 492];  % LC, PPN and nbM ROIs
selected_regions = [485, 486, 491, 492];  % LC and nbM ROIs
n_selected = length(selected_regions);

% Preallocate
mtd_selected = zeros(n_subjects, n_selected, n_ROIs, n_timepoints);

for i = 1:n_subjects
    for t = 1:n_timepoints
        A = matify(data(i, :, t)', n_ROIs);  % Reconstruct full 502x502 matrix
        mtd_selected(i, :, :, t) = A(selected_regions, :);  % Only keep rows from selected regions
    end

    % Print progress
    fprintf('Finished subject %d of %d\n', i, n_subjects);
end

mtd_lc_avg_csf4thvent = mean(mtd_selected(:,1:2,:,:),2); %take avg dFC of left & right LC
mtd_lc_avg_csf4thvent = reshape(mtd_lc_avg_csf4thvent,n_subjects,n_ROIs,n_timepoints);
mtd_nbm_avg_csf4thvent = mean(mtd_selected(:,3:4,:,:),2); %take avg dFC of left & right nbM
mtd_nbm_avg_csf4thvent = reshape(mtd_nbm_avg_csf4thvent,n_subjects,n_ROIs,n_timepoints);


mtd_lc_avg_del_csf4thvent = mtd_lc_avg_csf4thvent(delirium_sub_update,:,:);
mtd_lc_avg_health_csf4thvent = mtd_lc_avg_csf4thvent(health_sub_update,:,:);

% Avg. LC related MTD
for i=1:n_ROIs
   for k=1:n_timepoints
   [sig_mtd_lcavg_time_csf4thvent(i,k),pval_mtd_lcavg_time_csf4thvent(i,k)]=perm_code(mtd_lc_avg_del_csf4thvent(:,i,k),mtd_lc_avg_health_csf4thvent(:,i,k),1000);
   end
end

mean_lcavg_mtd_del = squeeze(mean(mtd_lc_avg_del_csf4thvent,1));
mean_lcavg_mtd_health = squeeze(mean(mtd_lc_avg_health_csf4thvent,1));
diff_mean_lcavg_mtd_csf4thvent = mean_lcavg_mtd_del - mean_lcavg_mtd_health; %any positive means greater coupling for Delirium

% Avg. nbm related MTD
mtd_nbm_avg_del_csf4thvent = mtd_nbm_avg_csf4thvent(delirium_sub_update,:,:);
mtd_nbm_avg_health_csf4thvent = mtd_nbm_avg_csf4thvent(health_sub_update,:,:);

for i=1:n_ROIs
   for k=1:n_timepoints
   [sig_mtd_nbmavg_time_csf4thvent(i,k),pval_mtd_nbmavg_time_csf4thvent(i,k)]=perm_code(mtd_nbm_avg_del_csf4thvent(:,i,k),mtd_nbm_avg_health_csf4thvent(:,i,k),1000);
   end
end
% FDR correction
pvals_flat = pval_mtd_nbmavg_time_csf4thvent(:);

mean_nbmavg_mtd_del = squeeze(mean(mtd_nbm_avg_del_csf4thvent,1));
mean_nbmavg_mtd_health = squeeze(mean(mtd_nbm_avg_health_csf4thvent,1));
diff_mean_nbmavg_mtd_csf4thvent = mean_nbmavg_mtd_del - mean_nbmavg_mtd_health;


% output for LDA Analysis -
% subject x time, ROIs
mtd_nbm_avg_reshape=permute(mtd_nbm_avg_csf4thvent,[1 3 2]);
mtd_nbm_avg_allsubstime = reshape(mtd_nbm_avg_reshape,117*210,502);
writematrix(mtd_nbm_avg_allsubstime,'mtd_nbm_avg_allsubstime.csv');

mtd_lc_avg_reshape=permute(mtd_lc_avg_csf4thvent,[1 3 2]);
mtd_lc_avg_allsubstime = reshape(mtd_lc_avg_reshape,117*210,502);
writematrix(mtd_lc_avg_allsubstime,'mtd_lc_avg_allsubstime.csv');

%% Compare between:

%original mtd_lc_avg - remove the subjects
mtd_lc_avg_orig_update=mtd_lc_avg;
mtd_lc_avg_orig_update(sub_remove_grp1,:,:)=[];

mtd_nbm_avg_orig_update=mtd_nbm_avg;
mtd_nbm_avg_orig_update(sub_remove_grp1,:,:)=[];

% t-test comparison across same subjects
% X, Y: subjects x regions x time
X = mtd_lc_avg_orig_update(:,1:400,:);
Y = mtd_lc_avg_csf4thvent(:,1:400,:);

[nSub, nReg, nTime] = size(X);

% Reshape to subjects x (regions*time)
X2 = reshape(X, nSub, []);
Y2 = reshape(Y, nSub, []);

% Paired t-test across subjects for each region-timepoint
[h, p] = ttest(X2, Y2);  % returns 1 x (nReg*nTime)

% Reshape back to region x time
p_rt = reshape(p, nReg, nTime);

% Optional: FDR correction across all tests
p_fdr = mafdr(p(:), 'BHFDR', true);
p_fdr_rt = reshape(p_fdr, nReg, nTime);
% effect size per sub
x = mean(mtd_lc_avg_orig_update(:,1:400,:), [2 3]); % 117x1
y = mean(mtd_lc_avg_csf4thvent(:,1:400,:), [2 3]);  % 117x1

Effect = meanEffectSize(x, y, Paired=true, Effect="cohen");
%----------------------------------------%
% Generate outputs of significance for plotting:
% just plot of the avg. significant over time (look at when they're different)
group_diff = diff_mean_lcavg_mtd_csf4thvent;
sig_matrix = sig_mtd_lcavg_time_csf4thvent;

group_diff = diff_mean_nbmavg_mtd_csf4thvent;
sig_matrix = sig_mtd_nbmavg_time_csf4thvent;
% group_diff is your original difference matrix [400 × timepoints]
masked_diff = group_diff .* sig_matrix;  % zero out non-significant entries

% Now average only across non-zero (i.e., significant) timepoints for each region
mean_sig_diff = sum(masked_diff, 2) ./ sum(sig_matrix, 2);  % [400 x 1]

% This handles division by zero (no significant timepoints) safely
mean_sig_diff(isnan(mean_sig_diff)) = 0;

% Assuming mean_sig_diff is [400 x 1]
filename = 'mean_sig_mtd_nbmavg.csv';

% Optional: Add region labels as row headers
% If you don't have labels, just write the data:
writematrix(mean_sig_diff, filename);

% Combine the significance and effect size threshold
sig_and_large = sig_matrix & (abs(group_diff) > 0.2);

% Mask the group differences
masked_diff = group_diff .* sig_and_large;

% Average across significant + large-effect timepoints only
mean_sig_diff = sum(masked_diff, 2) ./ sum(sig_and_large, 2);

% Handle divide-by-zero for regions with no significant large effects
mean_sig_diff(isnan(mean_sig_diff)) = 0;
writematrix(mean_sig_diff, 'mean_sig_mtd_nbmavg_thresholded.csv');

%just average across all timepoints - regardless
masked_diff = group_diff .* sig_matrix; 
mean_alltimepoint_diff = mean(masked_diff,2);
writematrix(mean_alltimepoint_diff, 'mean_sig_mtd_nbmavg_alltimepoints.csv');

% alternate threshold - threshold first & then take significant points
%set threshold:
threshold =0.2;
% Step 1: magnitude threshold (both + and -)
mag_mask = abs(group_diff) > threshold;

% Step 2: combine with significance mask
final_mask = mag_mask & sig_matrix;

% Step 3: apply mask
masked_diff = group_diff .* final_mask;

% Step 4: average across time (per ROI)
mean_sig_diff = sum(masked_diff, 2) ./ sum(final_mask, 2);
mean_sig_diff(isnan(mean_sig_diff)) = 0;

writematrix(mean_sig_diff, 'mean_sig_mtd_nbmavg_thresholdfirst.csv');


%----------------------------------------%
%% RECALCULATE FC
%----------------------------------------%
output_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';

files = dir(fullfile(output_dir, '*_csf4thvent_regress_timeseries.mat'));
if isempty(files)
    error('No timeseries files found in %s', output_dir);
end

% Determine ROI count from first file - avg of timeseries approach (not
% recommended)
tmp = load(fullfile(output_dir, files(1).name));
if isfield(tmp, 'cleaned_timeseries')
    ts0 = tmp.cleaned_timeseries;
elseif isfield(tmp, 'timeseries')
    ts0 = tmp.timeseries;
else
    error('No timeseries variable found in %s', files(1).name);
end

n_rois = size(484, 2);
n_targets = 2;
fc_aas_csf4thvent_regress = zeros(length(files), 484, 484);
sub_ids = strings(length(files), 1);

for ii = 1:length(files)
    subject_file = files(ii).name;
    split = strsplit(subject_file, '_');
    sub_ids(ii) = string(split{1});

    data = load(fullfile(output_dir, subject_file));
    if isfield(data, 'cleaned_timeseries')
        ts = data.cleaned_timeseries;
    elseif isfield(data, 'timeseries')
        ts = data.timeseries;
    else
        error('No timeseries variable found in %s', subject_file);
    end

    % Remove first/last 5 timepoints for noise
    if size(ts, 1) > 10
        ts = ts(6:end-5, :);
    end
    % average for aas - nbM & LC
    avg_ts_aas(:,1)=mean(ts(:,485:486),2); %LC
    avg_ts_aas(:,2)=mean(ts(:,491:492),2); %nBM

    %combine the timeseries together
    ts_no_aas = ts(:,1:482); % all non-aas
    ts_all = [ts_no_aas,avg_ts_aas]; % ts_all = 231 x 484

    % Functional connectivity (Pearson correlation)
    ts_corr = corr(ts_all);
    fc_aas_csf4thvent_regress(ii, :, :) = ts_corr;
    clear ts avg_ts_aas ts_no_aas ts_all
end


% remove relevant subjects
fc_aas_csf4thvent_regress(88,:,:)=[];
fc_aas_csf4thvent_regress(sub_remove_grp1,:,:)=[];

%group for delirium & health
fc_del_aas_csf4thvent = fc_aas_csf4thvent_regress(delirium_sub_update,:,:);
fc_health_aas_csf4thvent =fc_aas_csf4thvent_regress(health_sub_update,:,:);

%Grp differences - permuted
fc_aas_edges_del = reshape(fc_del_aas_csf4thvent,30,484*484);
fc_aas_edges_health = reshape(fc_health_aas_csf4thvent,87,484*484);
% run permutation across fc edges for aas - from AAS_analysis.m script (lines 573)
nEdges= length(fc_aas_edges_health);
for i=1:nEdges
    [sig_fc_edges1(i,:),pval_fc_edges1(i,:)]=perm_code(fc_aas_edges_del(:,i),fc_aas_edges_health(:,i),1000);
end
sig_fc_edges1 = reshape(sig_fc_edges1,484,484);
sig_fc_del_aas1 = squeeze(mean(fc_del_aas_csf4thvent,1)).*sig_fc_edges1;
sig_fc_health_aas1 = squeeze(mean(fc_health_aas_csf4thvent,1)).*sig_fc_edges1;
sig_diff_fc_del_health_aas1_csf4thvent = sig_fc_del_aas1 - sig_fc_health_aas1;
%sig diff for the 2 grps
diff_fc_del_health_lc_csf4thvent1 = diff_fc_del_health_aas1_csf4thvent(1:400,483);
writematrix(diff_fc_del_health_lc_csf4thvent1,'sig_diff_lc_fc_update1.csv')

diff_fc_del_health_nbm_csf4thvent1 = diff_fc_del_health_aas1_csf4thvent(1:400,484);
writematrix(diff_fc_del_health_nbm_csf4thvent1,'sig_diff_nbm_fc_update1.csv')

% extract LC/nbM to cortex FCs: (non-sig)
fc_all_lc_csf4thvent1 = fc_aas_csf4thvent_regress(:,1:400,483);
writematrix(fc_all_lc_csf4thvent1,'avg_LC_fc_400_all_subs_meants.csv');
fc_all_nbm_csf4thvent1 = fc_aas_csf4thvent_regress(:,1:400,484);
writematrix(fc_all_nbm_csf4thvent1,'avg_nbm_fc_400_all_subs_meants.csv');

% take the average for each - then difference than sig only
mean_fc_lc_csf4thvent1_del = squeeze(mean(fc_del_aas_csf4thvent(:,1:400,483)));
mean_fc_lc_csf4thvent1_health = squeeze(mean(fc_health_aas_csf4thvent(:,1:400,483)));

mean_fc_nbm_csf4thvent1_del = squeeze(mean(fc_del_aas_csf4thvent(:,1:400,484)));
mean_fc_nbm_csf4thvent1_health = squeeze(mean(fc_health_aas_csf4thvent(:,1:400,484)));

% take the difference
diff_fc_del_health_nbm_csf4thvent1 = mean_fc_nbm_csf4thvent1_del - mean_fc_nbm_csf4thvent1_health;
diff_fc_del_health_lc_csf4thvent1 = mean_fc_lc_csf4thvent1_del - mean_fc_lc_csf4thvent1_health;
% multiplication for sig
sig_diff_fc_del_health_nbm_csf4thvent1 = diff_fc_del_health_nbm_csf4thvent1.*sig_fc_edges1(1:400,484)';
writematrix(sig_diff_fc_del_health_nbm_csf4thvent1','sig_diff_nbm_fc_update2.csv')
sig_diff_fc_del_health_lc_csf4thvent1 = diff_fc_del_health_lc_csf4thvent1.*sig_fc_edges1(1:400,483)';
writematrix(sig_diff_fc_del_health_lc_csf4thvent1','sig_diff_lc_fc_update2.csv')


%% Alternate calculation
% Taking the flattened lower-triangle for significants

load()

% remove relevant subjects
fc_aas_csf4thvent_regress = fc_all_csf4thvent_regress;

fc_aas_csf4thvent_regress(88,:,:)=[];
fc_aas_csf4thvent_regress(sub_remove_grp1,:,:)=[];


nROIs = size(fc_aas_csf4thvent_regress,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_aas_csf4thvent_regress,117,nROIs*nROIs);
fc_edges_aas = flat_fc_all(:,template);
fc_edges_del = fc_edges_aas(delirium_sub_update,:); %fc edges lower tril delirium subs
fc_edges_health = fc_edges_aas(health_sub_update,:);
nEdges = size(fc_edges_aas,2);


fc_


fc_aas_csf4thvent_regress_only = fc_aas_csf4thvent_regress(:,[485:486,491:492],:);

fc_ass_csf4thvent_regress_only_nbm = squeeze(mean(fc_aas_csf4thvent_regress_only(:,[3:4],:),2));
fc_ass_csf4thvent_regress_only_lc = squeeze(mean(fc_aas_csf4thvent_regress_only(:,[1:2],:),2));

writematrix(fc_ass_csf4thvent_regress_only_lc,'avg_fc_lc_all_subs.csv');
writematrix(fc_ass_csf4thvent_regress_only_nbm,'avg_fc_nbm_all_subs.csv');

fc_edges_aas_nuclei_only = reshape(fc_aas_csf4thvent_regress_only,117,4*502);

fc_edges_aas_nuclei_del = fc_edges_aas_nuclei_only(delirium_sub_update,:);
fc_edges_aas_nuclei_health = fc_edges_aas_nuclei_only(health_sub_update,:);
nEdges=size(fc_edges_aas_nuclei_del,2);
for i=1:nEdges
    [sig_fc_edges_nuclei(i,:),pval_fc_edges_nuclei(i,:)]=perm_code(fc_edges_aas_nuclei_del(:,i),fc_edges_aas_nuclei_health(:,i),1000);
end

sig_fc_del_nuclei = sig_fc_edges_nuclei'.*squeeze(mean(fc_edges_aas_nuclei_del,1));
sig_fc_health_nuclei = sig_fc_edges_nuclei'.*squeeze(mean(fc_edges_aas_nuclei_health,1));

% calculate the deltas:
diff_fc_del_health_nuclei = sig_fc_del_nuclei - sig_fc_health_nuclei;

% separate the significant indices from there:
diff_fc_del_health_reshape = reshape(diff_fc_del_health_nuclei,4,502);


diff_fc_del_health_nosig = squeeze(mean(fc_edges_aas_nuclei_del,1)) - squeeze(mean(fc_edges_aas_nuclei_health,1));
diff_fc_del_health_sig = diff_fc_del_health_nosig.*sig_fc_edges_nuclei';

diff_fc_lc = diff_fc_del_health_reshape(1:2,1:400)';
diff_fc_nbm = diff_fc_del_health_reshape(3:4,1:400)';


% combine left & right together
A = diff_fc_lc(:,1);
B = diff_fc_lc(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
RB_surf_schaef_pink(data(1:400,1),'combine_lc_diff_sig_fc')
%min red and max green
RB_surf_schaef_color(data(1:400,1),'combine_lc_diff_sig_fc',[196 64 48],[12 125 121])
writematrix(data,'lc_combined_fc_sig_plot.csv')
%figure - ACh connectivity difference across brain

RB_surf_schaef_pink(diff_fc_nbm(1:400,1),'left_nbm_diff_sig_fc')
RB_surf_schaef_pink(diff_fc_nbm(1:400,2),'right_nbm_diff_sig_fc')
A = diff_fc_nbm(:,1);
B = diff_fc_nbm(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
RB_surf_schaef_color(data(1:400,1),'combine_nbm_diff_sig_fc',[196 64 48],[12 125 121])
writematrix(data,'nbm_combined_fc_sig_plot.csv')



% separate into nuclei:
diff_fc_del_health_aas2_csf4thvent_lc=diff_fc_del_health_aas_csf4thvent(1:400,483);
writematrix(diff_fc_del_health_aas_csf4thvent_lc,'sig_diff_lc_fc_csf4thvent_meantsnuclei.csv');

diff_fc_del_health_aas2_csf4thvent_nbm=diff_fc_del_health_aas_csf4thvent(1:400,484);
writematrix(diff_fc_del_health_aas_csf4thvent_nbm,'sig_diff_nbm_fc_csf4thvent_meantsnuclei.csv');

% significant difference in FC LC & nbM nuclei
sig_fc_edges = reshape(sig_fc_edges,484,484);

sig_fc_del_aas = squeeze(mean(fc_del_aas_csf4thvent,1)).*sig_fc_edges;
sig_fc_health_aas = squeeze(mean(fc_health_aas_csf4thvent,1)).*sig_fc_edges;

diff_fc_del_health_aas_csf4thvent = sig_fc_del_aas - sig_fc_health_aas;

% separate into nuclei:
diff_fc_del_health_aas_csf4thvent_lc=diff_fc_del_health_aas_csf4thvent(1:400,483);
writematrix(diff_fc_del_health_aas_csf4thvent_lc,'sig_diff_lc_fc_csf4thvent_meantsnuclei.csv');

diff_fc_del_health_aas_csf4thvent_nbm=diff_fc_del_health_aas_csf4thvent(1:400,484);
writematrix(diff_fc_del_health_aas_csf4thvent_nbm,'sig_diff_nbm_fc_csf4thvent_meantsnuclei.csv');



% load in the fc for the csf 4ht ventricle extraction
load('/Users/ntaylor/Desktop/AASDelirium_timeseries_csf4thvent/fc_all_csf4thvent_regress.mat')
% OLD METHOD - Extract FC values from the nuclei and then average them
% together

% remove the 1.5T subjects
%additional removal of 93 row - as this isn't calculated in the graph
%theory analysis
fc_all_csf4thvent_regress(88,:,:)=[]; %sub-164
fc_all_csf4thvent_regress(sub_remove_grp1,:,:)=[];

fc_all_csf4thvent_del = fc_all_csf4thvent_regress(delirium_sub_update,:,:);
fc_all_csf4thvent_health = fc_all_csf4thvent_regress(health_sub_update,:,:);

% flatten across lower triangle 
nROIs = size(fc_all_csf4thvent_regress,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_all_csf4thvent_regress,117,nROIs*nROIs);
fc_edges_all2 = flat_fc_all(:,template);
fc_edges_del2 = fc_edges_all2(delirium_sub_update,:); %fc edges lower tril delirium subs
fc_edges_health2 = fc_edges_all2(health_sub_update,:);
nEdges = size(fc_edges_all2,2);

% run permutation across fc edges - this is the correct way! - 9/09/24 -
% same way in AAS_analysis.m (line 30)
for i=1:nEdges
    [sig_fc_edges2(i,:),pval_fc_edges2(i,:)]=perm_code(fc_edges_del2(:,i),fc_edges_health2(:,i),1000);
end

sig_fc_mat = matify(sig_fc_edges2,nROIs); % rematify
mean_fc_all_csf4thvent_del = squeeze(mean(fc_all_csf4thvent_del));
mean_fc_all_csf4thvent_health = squeeze(mean(fc_all_csf4thvent_health));
sig_fc_all_csf4thvent_del2 = squeeze(mean(fc_all_csf4thvent_del)).*sig_fc_mat;
sig_fc_all_csf4thvent_health2 = squeeze(mean(fc_all_csf4thvent_health)).*sig_fc_mat;

diff_sig_fc_all_csf4thvent_del_health = sig_fc_all_csf4thvent_del2 - sig_fc_all_csf4thvent_health2;
% create difference first
diff_fc_all_csf4thvent_del_health = mean_fc_all_csf4thvent_del - mean_fc_all_csf4thvent_health;
% then significance ROIs only
diff_sig2_fc_all_csf4thvent_del_health = diff_fc_all_csf4thvent_del_health.*sig_fc_mat;


% create variables for plotting:
diff_fc_lc_csf4thvent=diff_sig_fc_all_csf4thvent_del_health(1:400,485:486);
% combine left & right together LC hemisphere differences:
A = diff_fc_lc_csf4thvent(:,1);
B = diff_fc_lc_csf4thvent(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
writematrix(data,'lc_combined_fc_sig_diff_csf4thvent.csv')
%nbM
diff_fc_nbM_csf4thvent=diff_sig_fc_all_csf4thvent_del_health(1:400,491:492);
% combine left & right together nbM hemisphere differences:
A = diff_fc_nbM_csf4thvent(:,1);
B = diff_fc_nbM_csf4thvent(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
writematrix(data,'nbM_combined_fc_sig_diff_csf4thvent.csv')



% generate significant just for nbM & LC
% LC indexed - 485:486
% nbM index - 491:492
% just LC/nbM

fc_delirium_lc_csf4thvent = fc_all_csf4thvent_del(:,:,485:486);
fc_health_lc_csf4thvent = fc_all_csf4thvent_health(:,:,485:486);
%nbm
fc_delirium_nbm_csf4thvent = fc_all_csf4thvent_del(:,:,491:492);
fc_health_nbm_csf4thvent = fc_all_csf4thvent_health(:,:,491:492);

fc_delirium_nbm_csf4thvent_flat = reshape(fc_delirium_nbm_csf4thvent,30,502*2);
fc_health_nbm_csf4thvent_flat = reshape(fc_health_nbm_csf4thvent,87,502*2);

nEdges=502*2;
for i=1:nEdges
    [sig_fc_aas_only(i,:),pval_fc_aas_only(i,:)]=perm_code(fc_delirium_nbm_csf4thvent_flat(:,i),fc_health_nbm_csf4thvent_flat(:,i),500);
end
sig_fc_aas_reshape = reshape(sig_fc_aas_only,502,2);

sig_fc_nbm_del = squeeze(mean(fc_delirium_nbm_csf4thvent,1)).*sig_fc_aas_reshape;
sig_fc_nbm_health = squeeze(mean(fc_health_nbm_csf4thvent,1)).*sig_fc_aas_reshape;

sig3_diff_fc_nbm_del_health = sig_fc_nbm_del - sig_fc_nbm_health;
% combine left & right together nbM hemisphere differences:
A = sig3_diff_fc_nbm_del_health(1:400,1);
B = sig3_diff_fc_nbm_del_health(1:400,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
writematrix(data,'nbM_combined_fc_sig3_diff_csf4thvent.csv')




writematrix(sig_diff_nbm_fc,'sig_diff_nbm_fc_whole.csv'); 
writematrix(sig_diff_lc_fc,'sig_diff_lc_fc_whole.csv'); 


% Fisher z-transformed FC calculation first & then averaged for the nuclei



% flatten across lower triangle 
nROIs = size(fc_all_csf4thvent_regress,3);
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
flat_fc_all = reshape(fc_all_csf4thvent_regress,117,nROIs*nROIs);
fc_edges_all2 = flat_fc_all(:,template);
% need to fisher z- transform data:
z_fc_edges_all2 = atanh(fc_edges_all2);

z_fc_edges_del2 = z_fc_edges_all2(delirium_sub_update,:); %fc edges lower tril delirium subs
z_fc_edges_health2 = z_fc_edges_all2(health_sub_update,:);
nEdges = size(z_fc_edges_all2,2);

% same way in AAS_analysis.m (line 30)
for i=1:nEdges
    [sig_fc_edges2_z(i,:),pval_fc_edges2_z(i,:)]=perm_code(z_fc_edges_del2(:,i),z_fc_edges_health2(:,i),1000);
end

sig_fc_mat_z = matify(sig_fc_edges2_z,nROIs); % rematify

mean_z_fc_csf4thvent_health = mean(z_fc_edges_health2);
mean_z_fc_csf4thvent_del = mean(z_fc_edges_del2);

diff_sig_fc_all_csf4thvent_del_health = sig_fc_all_csf4thvent_del2 - sig_fc_all_csf4thvent_health2;
% create difference first
diff_fc_all_csf4thvent_del_health = mean_fc_all_csf4thvent_del - mean_fc_all_csf4thvent_health;
% then significance ROIs only
diff_sig2_fc_all_csf4thvent_del_health = diff_fc_all_csf4thvent_del_health.*sig_fc_mat;


% create variables for plotting:
diff_fc_lc_csf4thvent=diff_sig_fc_all_csf4thvent_del_health(1:400,485:486);
% combine left & right together LC hemisphere differences:
A = diff_fc_lc_csf4thvent(:,1);
B = diff_fc_lc_csf4thvent(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
writematrix(data,'lc_combined_fc_sig_diff_csf4thvent.csv')
%nbM
diff_fc_nbM_csf4thvent=diff_sig_fc_all_csf4thvent_del_health(1:400,491:492);
% combine left & right together nbM hemisphere differences:
A = diff_fc_nbM_csf4thvent(:,1);
B = diff_fc_nbM_csf4thvent(:,2);
data = (A == 0) .* B + (B == 0) .* A + (A ~= 0 & B ~= 0) .* ((A + B) / 2);
writematrix(data,'nbM_combined_fc_sig_diff_csf4thvent.csv')


%----------------------------------------%
%% AAS Cross-Corr
%----------------------------------------%

%generate a new script for this recalculation - "AAS_Cross_Corr_update.m"


