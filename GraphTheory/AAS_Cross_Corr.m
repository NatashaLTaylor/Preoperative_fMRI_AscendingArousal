%% Cross-Correlation of AAS Nuclei to Whole-Brain Connectivity

load('timeseries_schaef_all_preop.csv');

% Load time-series data (T time points x N brain regions)
% Assume data is stored in a matrix 'data' of size T x N

cd '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400'
filename = dir("*timeseries.mat");

max_timepoints = 221;
num_rois = 502;  % number of ROIs
ts_all = nan(length(filename), num_rois, max_timepoints);  % preallocate with NaNs
sub_ids = cell(length(filename), 1);  % store subject IDs

for ii = 1:length(filename)
    % Get 'sub-...' from filename
    subject_file = filename(ii).name;
    split = strsplit(subject_file, '_');
    subnum = split{1};
    sub_ids{ii} = subnum;

    % Load and trim ts
    load(subject_file);  % assumes variable 'ts' exists
    ts = ts(6:end-5, :);  % remove first and last 5 timepoints

    % Check length
    if size(ts, 1) == max_timepoints
        ts_all(ii, :, :) = ts';  % transpose to roi x time
    else
        ts_all(ii, :, :) = NaN;  % entire time series replaced with NaNs
    end
end
%now remove relevant IDs
sub_remove_clinical = [42,43,59,63,74,93,103,104];
ts_all(sub_remove_clinical,:,:)=[];



%% Running cross-correlation in just the fc?
%Define seed region index (change as needed)
selected_regions = [485, 486, 489, 490, 491, 492];  % LC, PPN and nbM ROIs
n_selected = length(selected_regions);

seed_idx = 485; % Example: first region as seed

% Define range of time lags
max_lag = 10; % Maximum lag to consider (e.g., +/- 10 time points)
lags = -max_lag:max_lag;

% Get seed region time series
seed_ts = data(:, seed_idx);

% Initialize cross-correlation matrix
num_regions = size(data, 2);
cross_corr = zeros(length(lags), num_regions);

% Compute cross-correlation for each region - just for BOLD time-series
for region_idx = 1:num_regions
    if region_idx == seed_idx
        continue; % Skip self-correlation
    end
    region_ts = data(:, region_idx);
    
    % Compute cross-correlation
    [xc, lags_out] = xcorr(seed_ts, region_ts, max_lag, 'coeff');
    
    % Store only the relevant range
    cross_corr(:, region_idx) = xc;
end


% Plot cross-correlation
figure;
imagesc(lags, 1:num_regions, cross_corr');
colorbar;
ylabel('Brain Regions');
xlabel('Time Lag');
title('Cross-Correlation of Seed Region with Other Regions');


% can I look at the cross correlation between BOLD of LC and nbM & PC
% integration segregation?

%% Change in BOLD LC, nbM & PPN vs PC values
% load in relevant PC values
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/participation_coeff_all.mat')
% part_all = 125 x 502 x 210

% ===============================
%% 2. Mean PC (over all ROIs) related to changes in difference in LC, PPN, nbm
% ===============================


% turn into a structure
seed_names = {'LC','PPN','nbM','ACh'};
results_rawts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [],...
    'std_ts', [], ...
    'xcf', [], ...
    'lags',[], ...
    'max_corr',[], ...
    'max_lag',[]);

window_size = 5;  % Define your window size

for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,7:end-5);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    ach_ts = mean(ts_compare([489:492],:), 1);

    seed_ts_list = {lc_ts, ppn_ts, nbm_ts, ach_ts};
    % Loop through each seed
    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};    
            [xcf, lags] = crosscorr(mean_pc, seed_ts, NumLags=window_size);
            [~, I] = max(abs(xcf));
            maxLags = I;
            maxCorrs = xcf(I);
        % Save all results in the structure
        results_rawts(sub, s).subject   = sub;
        results_rawts(sub, s).seed_name = seed_names{s};
        results_rawts(sub, s).timeseries = seed_ts;
        results_rawts(sub, s).std_ts    = std_ts;
        results_rawts(sub, s).xcf  = xcf;
        results_rawts(sub, s).lags      = lags;
        results_rawts(sub, s).max_corr  = maxCorr;
        results_rawts(sub, s).max_lag   = maxLag;
    end
        fprintf("Subject finished %d\n", sub);
end

% plot all 4 - with colour grouping
nSubs = size(results_rawts, 1);    % e.g., 119
nSeeds = size(results_rawts, 2);   % e.g., 4

%% significance testing
% look at separating by the group 
% Initialize the structure before the loop
xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_diffrawts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_del(i, :) = this_entry.xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_health(i, :) = this_entry.xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Store in structure
    xcf_by_group(seed_idx).seed_name = seed_names{seed_idx};
    xcf_by_group(seed_idx).lags = lags;

    xcf_by_group(seed_idx).delirium = all_xcf_del;
    xcf_by_group(seed_idx).health = all_xcf_health;

    xcf_by_group(seed_idx).mean_delirium = group_mean_del;
    xcf_by_group(seed_idx).sem_delirium  = group_sem_del;
    xcf_by_group(seed_idx).mean_health  = group_mean_health;
    xcf_by_group(seed_idx).sem_health   = group_sem_health;
end


% filter through - normality test and then relevant significance test run

for seed_idx = 1:nSeeds
    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);
    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);
    test_type = strings(1, nLags);  % Track which test was used

    for lag_idx = 1:nLags
        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        if length(group1) > 1 && length(group2) > 1
            % Normality check
            is_normal1 = lillietest(group1) == 0;
            is_normal2 = lillietest(group2) == 0;

            if is_normal1 && is_normal2
                % Use t-test
                [~, p, ~, stats] = ttest2(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = stats.tstat;
                test_type(lag_idx) = "ttest";
            else
                % Use non-parametric Mann–Whitney U test
                p = ranksum(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = NaN;  % t-stat not applicable
                test_type(lag_idx) = "ranksum";
            end
        end
    end

    % Store in structure
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;
    xcf_by_group(seed_idx).test_type = test_type;
end


% run permuted p-val
% sub x lag (randomise grp IDs to determine significance)
sig_permute= nan(nSeeds,length(lags));
pval_permute= nan(nSeeds,length(lags));
for seed_idx = 1:nSeeds
    del_grp_crosscor= xcf_by_group(seed_idx).delirium;
    health_grp_crosscor= xcf_by_group(seed_idx).health;
    for i=1:length(lags)
        [sig_permute(seed_idx,i),pval_permute(seed_idx,i)]=perm_code(del_grp_crosscor(:,i),health_grp_crosscor(:,i),1000);
    end
    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;
end
% check if any are significant
for seed_idx = 1:nSeeds
    p_vals = xcf_by_group(seed_idx).p_vals;
    lags = xcf_by_group(seed_idx).lags;
    seed_name = xcf_by_group(seed_idx).seed_name;
    
    sig_lags = find(p_vals < 0.05);

    if ~isempty(sig_lags)
        fprintf("Seed %s: %d significant lag(s) found (p < 0.05)\n", ...
            seed_name, length(sig_lags));
        fprintf("  Lags: %s\n", mat2str(lags(sig_lags)));
        fprintf("  p-values: %s\n", mat2str(p_vals(sig_lags), 3));  % rounded to 3 decimals
    else
        fprintf("Seed %s: No significant lags (p < 0.05)\n", seed_name);
    end
end


%rename so variables don't get confused
xcf_rawts_grp_stats = xcf_by_group;


%% figures
figure;
set(gcf, 'color', 'w');  % White background

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.4039, 0.6392, 0.6235;      % Dark green for seed 2 (PPN)
    0.2627, 0.4196, 0.3725;      % Medium green for seed 3 (nbM)
    0.1176, 0.5098, 0.4941       % Light green for seed 4 (ACh)
];

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_rawts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_rawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_mean_xcf(sub, :) = this_entry.xcf;
            plot(lags, this_entry.xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_rawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

filename = fullfile(pwd, 'crosscorr_rawts_allsub_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');


% just plot the mean for simplicity
nSubs = size(results_rawts, 1);    % e.g., 119
nSeeds = size(results_rawts, 2);   % e.g., 4


figure;
set(gcf, 'color', 'w');  % White background

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.4039, 0.6392, 0.6235;      % Dark green for seed 2 (PPN)
    0.2627, 0.4196, 0.3725;      % Medium green for seed 3 (nbM)
    0.1176, 0.5098, 0.4941       % Light green for seed 4 (ACh)
];

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_rawts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_rawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_mean_xcf(sub, :) = this_entry.xcf;
        end
    end

    % Plot group mean XCF
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    plot(lags, group_mean, 'Color', seed_colors(seed_idx, :), 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_rawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean XCF');
end


filename = fullfile(pwd, 'crosscorr_rawts_mean_xcf.svg');
print(gcf, filename, '-dsvg');


%plotting the group means overlayed
figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_rawts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_rawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_del(i, :) = this_entry.xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_rawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_health(i, :) = this_entry.xcf;
        end
    end

    % Compute group means & sem
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Plot shaded error region for health
    fill([lags, fliplr(lags)], ...
         [group_mean_health + group_sem_health, fliplr(group_mean_health - group_sem_health)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % Plot shaded error region for delirium
    fill([lags, fliplr(lags)], ...
         [group_mean_del + group_sem_del, fliplr(group_mean_del - group_sem_del)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot both lines with same color but different line styles
    h1 = plot(lags, group_mean_del, '--', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % dashed
    h2 = plot(lags, group_mean_health, '-', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % solid

    % Titles and labels
    title(strrep(results_rawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean Cross-corr');

    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

% Save figure as .svg
filename = fullfile(pwd, 'crosscorr_rawts_mean_xcf_del_vs_health_with_sem.svg');
print(gcf, filename, '-dsvg');





% ===============================
%% 2. Lops through each diff in BOLD timeseries and relate to
% ===============================
% avg. PC for the time point
seed_names = {'LC -PPN','LC - nbM', 'LC - ACh', 'PPN - LC', 'nbM - LC', 'ACh - LC', 'LC + nbM'};
results_diffrawts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [],...
    'std_ts', [],...
    'xcf', [], ...
    'lags',[], ...
    'max_corr',[], ...
    'max_lag',[]);

window_size = 5;  % Define your window size

for sub = 1:nSubs
    % Get time series for this subject
    ts = squeeze(ts_all(sub, :, :));            % [time x region]
    ts_compare = ts(:,7:end-5);                % trim edges
    part = squeeze(part_all(sub,:,:));
    mean_pc = mean(part,1);
    % Define seed time series (adjusting for trimmed indices)
    lc_ts  = mean(ts_compare([485,486],:), 1);
    ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    ach_ts = mean(ts_compare([489:492],:),1);
    % take differences:
    diff_lc_ppn = lc_ts - ppn_ts;
    diff_lc_nbm = lc_ts - nbm_ts;
    diff_lc_ach = lc_ts - ach_ts;
    diff_ppn_lc = ppn_ts - lc_ts;
    diff_nbm_lc = nbm_ts - lc_ts;
    diff_ach_lc = ach_ts - lc_ts;
    lcnbm_combine = lc_ts + nbm_ts;

    seed_ts_list = {diff_lc_ppn,diff_lc_nbm, diff_lc_ach, diff_ppn_lc, diff_nbm_lc, diff_ach_lc,lcnbm_combine};
    % Loop through each seed
    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};    
        std_ts = std(seed_ts);
            [xcf, lags] = crosscorr(mean_pc, seed_ts, NumLags=window_size);
            [~, I] = max(abs(xcf));
            maxLag = I;
            maxCorr = xcf(I);
        % Save all results in the structure
        results_diffrawts(sub, s).subject   = sub;
        results_diffrawts(sub, s).seed_name = seed_names{s};
        results_diffrawts(sub, s).timeseries = seed_ts;
        results_diffrawts(sub, s).std_ts = std_ts;
        results_diffrawts(sub, s).xcf  = xcf;
        results_diffrawts(sub, s).lags  = lags;
        results_diffrawts(sub, s).max_corr = maxCorr;
        results_diffrawts(sub, s).max_lag  = maxLag;
    end
        fprintf("Subject finished %d\n", sub);
end



% is there a significant difference between the groups

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_diffrawts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_del(i, :) = this_entry.xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_health(i, :) = this_entry.xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

end


% look at separating by the group 
% Initialize the structure before the loop
xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_diffrawts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_del(i, :) = this_entry.xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_health(i, :) = this_entry.xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Store in structure
    xcf_by_group(seed_idx).seed_name = seed_names{seed_idx};
    xcf_by_group(seed_idx).lags = lags;

    xcf_by_group(seed_idx).delirium = all_xcf_del;
    xcf_by_group(seed_idx).health = all_xcf_health;

    xcf_by_group(seed_idx).mean_delirium = group_mean_del;
    xcf_by_group(seed_idx).sem_delirium  = group_sem_del;
    xcf_by_group(seed_idx).mean_health  = group_mean_health;
    xcf_by_group(seed_idx).sem_health   = group_sem_health;
end


% filter through - normality test and then relevant significance test run

for seed_idx = 1:nSeeds
    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);
    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);
    test_type = strings(1, nLags);  % Track which test was used

    for lag_idx = 1:nLags
        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        if length(group1) > 1 && length(group2) > 1
            % Normality check
            is_normal1 = lillietest(group1) == 0;
            is_normal2 = lillietest(group2) == 0;

            if is_normal1 && is_normal2
                % Use t-test
                [~, p, ~, stats] = ttest2(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = stats.tstat;
                test_type(lag_idx) = "ttest";
            else
                % Use non-parametric Mann–Whitney U test
                p = ranksum(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = NaN;  % t-stat not applicable
                test_type(lag_idx) = "ranksum";
            end
        end
    end

    % Store in structure
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;
    xcf_by_group(seed_idx).test_type = test_type;
end


% run permuted p-val
% sub x lag (randomise grp IDs to determine significance)
sig_permute= nan(nSeeds,length(lags));
pval_permute= nan(nSeeds,length(lags));
for seed_idx = 1:nSeeds
    del_grp_crosscor= xcf_by_group(seed_idx).delirium;
    health_grp_crosscor= xcf_by_group(seed_idx).health;
    for i=1:length(lags)
        [sig_permute(seed_idx,i),pval_permute(seed_idx,i)]=perm_code(del_grp_crosscor(:,i),health_grp_crosscor(:,i),1000);
    end
    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;
end
%rename so variables don't get confused
xcf_diffrawts_grp_stats = xcf_by_group;

% figures
nSubs = size(results_diffrawts, 1);    % e.g., 119
nSeeds = size(results_diffrawts, 2);   % e.g., 4

figure;
set(gcf, 'color', 'w');             % White background
% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.7608, 0.2392, 0.1804;
    0.7608, 0.2392, 0.1804;
    0.403921568627451, 0.6392156862745098, 0.6235294117647059;      % Dark green for seed 2 (PPN)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
    0.11764705882352941, 0.5098039215686274, 0.49411764705882355;% Light green for seed 4 (ACh)
    0.403921568627451, 0.49411764705882355,0.6470588235294118 %pale blue
];
for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject
    lags = results_diffrawts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_mean_xcf(sub, :) = this_entry.xcf;
            plot(lags, this_entry.xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_diffrawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Add a single legend to the last subplot (optional)
subplot(2, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');

filename = fullfile(pwd, 'crosscorr_diffrawts_allsubs_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');




% just plot the mean for simplicity
nSubs = size(results_diffrawts, 1);    % e.g., 119
nSeeds = size(results_diffrawts, 2);   % e.g., 4

figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_diffrawts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_mean_xcf(sub, :) = this_entry.xcf;
        end
    end

    % Plot group mean XCF
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    plot(lags, group_mean, 'Color', seed_colors(seed_idx, :), 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_diffrawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean XCF');
end

% Define filename and full path (optional)
filename = fullfile(pwd, 'crosscorr_diffrawts_mean_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');



% plot based upon grouping
group_subs = health_sub;
nSubs_group = length(group_subs);  % Number of subjects in group
nSeeds = size(results_diffrawts, 2);

figure;
set(gcf, 'color', 'w');  % White background

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;  % Red for seed 1 (LC)
    0.7608, 0.2392, 0.1804;
    0.7608, 0.2392, 0.1804;
    0.4039, 0.6392, 0.6235;  % Dark green for seed 2 (PPN)
    0.2627, 0.4196, 0.3725;  % Medium green for seed 3 (nbM)
    0.1176, 0.5098, 0.4941;  % Light green for seed 4 (ACh)
    0.4039, 0.4941, 0.6471;  % Pale blue
];

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject in group (to get lag range)
    lags = results_diffrawts(group_subs(1), seed_idx).lags;
    all_mean_xcf = nan(nSubs_group, length(lags));  % Preallocate

    for i = 1:nSubs_group
        sub = group_subs(i);  % Actual subject index
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_mean_xcf(i, :) = this_entry.xcf;
            plot(lags, this_entry.xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_diffrawts(group_subs(1), seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Add a single legend to the last subplot (optional)
subplot(2, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');

%plotting the group means overlayed
figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_diffrawts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_del(i, :) = this_entry.xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_diffrawts(sub, seed_idx);
        if ~isempty(this_entry.xcf)
            all_xcf_health(i, :) = this_entry.xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Plot shaded error region for health
    fill([lags, fliplr(lags)], ...
         [group_mean_health + group_sem_health, fliplr(group_mean_health - group_sem_health)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot shaded error region for delirium
    fill([lags, fliplr(lags)], ...
         [group_mean_del + group_sem_del, fliplr(group_mean_del - group_sem_del)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot both lines with same color but different line styles
    h1 = plot(lags, group_mean_del, '--', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % dashed for delirium
    h2 = plot(lags, group_mean_health, '-', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % solid for healthy

    % Titles and labels
    title(strrep(results_diffrawts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean Cross-corr');

    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

% Save figure as .svg
filename = fullfile(pwd, 'crosscorr_diffrawts_mean_xcf_del_vs_health_with_sem.svg');
print(gcf, filename, '-dsvg');








% ===============================
%% 3. just looking at the peaks in time-series -time-points with the peaks time-series, defined by 2* standard deviation (thresh)
% ===============================

seed_names = {'LC','PPN','nbM','ACh'};
% Preallocate structure array
results_peakts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [],...
    'peak_inds', [], ...
    'total_peaks',[],...
    'std_ts', [], ...
    'xcf', [], ...
    'mean_xcf',[],...
    'lags',[], ...
    'max_corr',[], ...
    'max_lag',[]);

window_size = 5;  % Define your window size

xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_health(i, :) = this_entry.mean_xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Store in structure
    xcf_by_group(seed_idx).seed_name = seed_names{seed_idx};
    xcf_by_group(seed_idx).lags = lags;

    xcf_by_group(seed_idx).delirium = all_xcf_del;
    xcf_by_group(seed_idx).health = all_xcf_health;

    xcf_by_group(seed_idx).mean_delirium = group_mean_del;
    xcf_by_group(seed_idx).sem_delirium  = group_sem_del;
    xcf_by_group(seed_idx).mean_health  = group_mean_health;
    xcf_by_group(seed_idx).sem_health   = group_sem_health;
end


% filter through - normality test and then relevant significance test run
%% significant testing
for seed_idx = 1:nSeeds
    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);
    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);
    test_type = strings(1, nLags);  % Track which test was used

    for lag_idx = 1:nLags
        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        if length(group1) > 1 && length(group2) > 1
            % Normality check
            is_normal1 = lillietest(group1) == 0;
            is_normal2 = lillietest(group2) == 0;

            if is_normal1 && is_normal2
                % Use t-test
                [~, p, ~, stats] = ttest2(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = stats.tstat;
                test_type(lag_idx) = "ttest";
            else
                % Use non-parametric Mann–Whitney U test
                p = ranksum(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = NaN;  % t-stat not applicable
                test_type(lag_idx) = "ranksum";
            end
        end
    end

    % Store in structure
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;
    xcf_by_group(seed_idx).test_type = test_type;
end


for seed_idx = 1:nSeeds
    p_vals = xcf_by_group(seed_idx).p_vals;
    lags = xcf_by_group(seed_idx).lags;
    seed_name = xcf_by_group(seed_idx).seed_name;
    
    sig_lags = find(p_vals < 0.05);

    if ~isempty(sig_lags)
        fprintf("Seed %s: %d significant lag(s) found (p < 0.05)\n", ...
            seed_name, length(sig_lags));
        fprintf("  Lags: %s\n", mat2str(lags(sig_lags)));
        fprintf("  p-values: %s\n", mat2str(p_vals(sig_lags), 3));  % rounded to 3 decimals
    else
        fprintf("Seed %s: No significant lags (p < 0.05)\n", seed_name);
    end
end

%% run permuted p-val
% sub x lag (randomise grp IDs to determine significance)
sig_permute= nan(nSeeds,length(lags));
pval_permute= nan(nSeeds,length(lags));
for seed_idx = 1:nSeeds
    del_grp_crosscor= xcf_by_group(seed_idx).delirium;
    health_grp_crosscor= xcf_by_group(seed_idx).health;
    for i=1:length(lags)
        [sig_permute(seed_idx,i),pval_permute(seed_idx,i)]=perm_code(del_grp_crosscor(:,i),health_grp_crosscor(:,i),1000);
    end
    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;
end
%rename so variables don't get confused
xcf_peakts_grp_stats = xcf_by_group;



%% Figures - Peak Timeseries (step 3)
figure;
set(gcf, 'Color', 'w');
seed_names = {'LC', 'PPN', 'nbM', 'ACh'};

for s = 1:4
    subplot(1, 4, s);
    hold on;

    % Get structure for this subject and seed
    this_entry = results_peakts(sub, s);
    
    % Extract data
    xcf_all = this_entry.xcf;      % [nPeaks x lags]
    lags    = this_entry.lags;

    % Plot each peak's XCF
    if ~isempty(xcf_all)
        for i = 1:size(xcf_all, 1)
            plot(lags, xcf_all(i, :));  % gray lines
        end

        % Plot average across peaks
        plot(lags, this_entry.mean_xcf, 'r', 'LineWidth', 2);  % mean in red
    end

    title(seed_names{s});
    xlabel('Lag');
    ylabel('Cross-corr');
    xlim([min(lags), max(lags)]);
    grid on;
end


% plot all 4 - with colour grouping
nSubs = size(results_peakts, 1);    % e.g., 119
nSeeds = size(results_peakts, 2);   % e.g., 4

figure;
set(gcf, 'color', 'w');             % White background
% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.403921568627451, 0.6392156862745098, 0.6235294117647059;      % Dark green for seed 2 (PPN)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
    0.11764705882352941, 0.5098039215686274, 0.49411764705882355       % Light green for seed 4 (ACh)
];
for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject
    lags = results_peakts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(sub, :) = this_entry.mean_xcf;
            plot(lags, this_entry.mean_xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Add a single legend to the last subplot (optional)
subplot(1, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');

filename = fullfile(pwd, 'crosscorr_peakts_allsub_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');

% just plot the mean for simplicity

figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(sub, :) = this_entry.mean_xcf;
        end
    end

    % Plot group mean XCF
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    plot(lags, group_mean, 'Color', seed_colors(seed_idx, :), 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean XCF');
end

filename = fullfile(pwd, 'crosscorr_peakts_mean_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');


%plotting the group means overlayed
figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_health(i, :) = this_entry.mean_xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));
        % Plot shaded error region for health
    fill([lags, fliplr(lags)], ...
         [group_mean_health + group_sem_health, fliplr(group_mean_health - group_sem_health)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot shaded error region for delirium
    fill([lags, fliplr(lags)], ...
         [group_mean_del + group_sem_del, fliplr(group_mean_del - group_sem_del)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot both lines with same color but different line styles
    h1 = plot(lags, group_mean_del, '--', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % dashed
    h2 = plot(lags, group_mean_health, '-', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % solid

    % Titles and labels
    title(strrep(results_peakts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean Cross-corr');

    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

% Save figure as .svg
filename = fullfile(pwd, 'crosscorr_peakts_mean_xcf_del_vs_health_wtih_sem.svg');
print(gcf, filename, '-dsvg');


% plot for delirium and non-delirium separate
% Define which subjects to include (e.g., based on clinical group)
group_subs = delirium_sub;  % <-- replace this with your actual subject IDs
group_subs = health_sub;
figure;
set(gcf, 'color', 'w');             % White background

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.4039, 0.6392, 0.6235;      % Dark green for seed 2 (PPN)
    0.2627, 0.4196, 0.3725;      % Medium green for seed 3 (nbM)
    0.1176, 0.5098, 0.4941       % Light green for seed 4 (ACh)
];

for seed_idx = 1:nSeeds
    subplot(1, 4, seed_idx);
    hold on;

    % Extract lags from the first subject
    lags = results_peakts(1, seed_idx).lags;
    all_mean_xcf = nan(length(group_subs), length(lags));  % Preallocate

    for i = 1:length(group_subs)
        sub = group_subs(i);
        this_entry = results_peakts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(i, :) = this_entry.mean_xcf;
            plot(lags, this_entry.mean_xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Optional: Add legend
subplot(1, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');




% look at total number of peaks for each time-series
nSubs = size(results_peakts, 1);
nSeeds = size(results_peakts, 2);

total_peaks_matrix = zeros(nSubs, nSeeds);

for sub = 1:nSubs
    for seed = 1:nSeeds
        total_peaks_matrix(sub, seed) = results_peakts(sub, seed).total_peaks;
    end
end

total_peaks_del = sum(total_peaks_matrix(delirium_sub, :), 2);
total_peaks_nondel = sum(total_peaks_matrix(health_sub, :), 2);

mean_peaks_del     = mean(total_peaks_matrix(delirium_sub, :), 1);
mean_peaks_nondel  = mean(total_peaks_matrix(health_sub, :), 1);
% Define colors
nondel_color = [253, 205, 154] / 255;  % #fdcd9a
del_color    = [174, 216, 230] / 255;  % #aed8e6

% Compute group stats
group_means = [
    mean(total_peaks_matrix(health_sub, :), 1);
    mean(total_peaks_matrix(delirium_sub, :), 1)
];  % [2 x nSeeds]

group_stds = [
    std(total_peaks_matrix(health_sub, :), 0, 1);
    std(total_peaks_matrix(delirium_sub, :), 0, 1)
];  % [2 x nSeeds]

% Create bar plot
figure;
set(gcf,'color','w')
hold on;
b = bar(group_means', 'grouped');
b(1).FaceColor = nondel_color;
b(2).FaceColor = del_color;

% Calculate x positions for error bars
num_seeds = length(seed_names);
group_width = min(0.8, 2 / 3);
x = zeros(2, num_seeds);
for i = 1:2
    x(i,:) = (1:num_seeds) - group_width/2 + (2*i-1) * group_width / 4;
end

% Add error bars
errorbar(x(1,:), group_means(1,:), group_stds(1,:), 'k.', 'LineWidth', 1);
errorbar(x(2,:), group_means(2,:), group_stds(2,:), 'k.', 'LineWidth', 1);

% Aesthetics
set(gca, 'XTick', 1:num_seeds, 'XTickLabel', seed_names, 'XTickLabelRotation', 45);
ylabel('Average Number of Peaks');
legend({'Non-Delirious', 'Delirious'}, 'Location', 'northeast');
title('Average Peaks per Seed by Group');
grid on;
hold off;

% significant
% t-tests:
p_vals = zeros(1, length(seed_names));
t_stats = zeros(1, length(seed_names));

for s = 1:length(seed_names)
    grp1 = total_peaks_matrix(health_sub, s);
    grp2 = total_peaks_matrix(delirium_sub, s);
    
    [~, p, ~, stats] = ttest2(grp1, grp2);
    
    p_vals(s) = p;
    t_stats(s) = stats.tstat;
end

% Display results
T = table(seed_names(:), t_stats(:), p_vals(:), 'VariableNames', {'Seed', 'T_stat', 'P_value'})
% Mann-Whitney Test - non-parametric assumption
p_vals_ranksum = zeros(1, length(seed_names));

for s = 1:length(seed_names)
    grp1 = total_peaks_matrix(health_sub, s);
    grp2 = total_peaks_matrix(delirium_sub, s);
    
    p = ranksum(grp1, grp2);  % Mann-Whitney U test
    p_vals_ranksum(s) = p;
end

% Display results
T_nonparam = table(seed_names(:), p_vals_ranksum(:), 'VariableNames', {'Seed', 'P_value_ranksum'})


% ===============================
%% 4. Look at peaks in delta time-series
% ===============================
seed_names = {'LC-PPN', 'LC-nbM', 'LC-ACh', 'PPN-LC', 'nbM-LC', 'ACh-LC', 'LC+nbM'};
% Preallocate structure array
results_peakdiffts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [],...
    'peak_inds', [], ...
    'total_peaks',[],...
    'std_ts', [], ...
    'xcf', [], ...
    'mean_xcf',[],...
    'lags',[], ...
    'max_corr',[], ...
    'max_lag',[]);

window_size = 5;  % Define your window size
%calculating across each sub
for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,7:end-5);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    ach_ts = mean(ts_compare([489:492],:), 1);

    lcppn_ts = lc_ts - ppn_ts;
    lcnbm_ts = lc_ts - nbm_ts;
    lcach_ts = lc_ts - ach_ts;
    ppnlc_ts = ppn_ts - lc_ts;
    nbmlc_ts = nbm_ts - lc_ts;
    achlc_ts = ach_ts - lc_ts;
    lcnbm_combine = lc_ts + nbm_ts;

    seed_ts_list = {lcppn_ts, lcnbm_ts, lcach_ts, ppnlc_ts, nbmlc_ts, achlc_ts, lcnbm_combine};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % Detect peaks above 2 SDs from mean
        stdThresh =2;
        [~,peak_inds] =  findpeaks(seed_ts,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(seed_ts)); %find peaks of acc
        %total peaks
        sum_pks = length(peak_inds);
        xcf_accum = [];

        for i = 1:length(peak_inds)
            idx = peak_inds(i);
            start_idx = idx - window_size;
            end_idx = idx + window_size;

            if start_idx < 1 || end_idx > length(seed_ts)
                continue;
            end

            seg_seed = seed_ts(start_idx:end_idx);
            seg_pc   = mean_pc(start_idx:end_idx);

            [xcf, lags] = crosscorr(seg_pc, seg_seed, NumLags=window_size);
            xcf_accum = [xcf_accum; xcf];
        end

        if ~isempty(xcf_accum)
            mean_xcf = mean(xcf_accum, 1);
            [~, I] = max(abs(mean_xcf));
            maxCorr = mean_xcf(I);
            maxLag = lags(I);
        else
            mean_xcf = nan(1, 2*window_size+1);  % assuming default lags from -10:10
            lags = -window_size:window_size;
            maxCorr = nan;
            maxLag = nan;
        end

        % Save all results in the structure
        results_peakdiffts(sub, s).subject   = sub;
        results_peakdiffts(sub, s).seed_name = seed_names{s};
        results_peakdiffts(sub, s).timeseries = seed_ts;
        results_peakdiffts(sub, s).peak_inds = peak_inds;
        results_peakdiffts(sub, s).total_peaks = sum_pks;
        results_peakdiffts(sub, s).std_ts    = std_ts;
        results_peakdiffts(sub, s).xcf  = xcf_accum;
        results_peakdiffts(sub, s).mean_xcf  = mean_xcf;
        results_peakdiffts(sub, s).lags      = lags;
        results_peakdiffts(sub, s).max_corr  = maxCorr;
        results_peakdiffts(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

% is there a significant difference between the groups

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_health(i, :) = this_entry.mean_xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));
end


% look at separating by the group 
% Initialize the structure before the loop
xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_health(i, :) = this_entry.mean_xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

    % Store in structure
    xcf_by_group(seed_idx).seed_name = seed_names{seed_idx};
    xcf_by_group(seed_idx).lags = lags;

    xcf_by_group(seed_idx).delirium = all_xcf_del;
    xcf_by_group(seed_idx).health = all_xcf_health;

    xcf_by_group(seed_idx).mean_delirium = group_mean_del;
    xcf_by_group(seed_idx).sem_delirium  = group_sem_del;
    xcf_by_group(seed_idx).mean_health  = group_mean_health;
    xcf_by_group(seed_idx).sem_health   = group_sem_health;
end


% filter through - normality test and then relevant significance test run
%% significant testing
for seed_idx = 1:nSeeds
    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);
    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);
    test_type = strings(1, nLags);  % Track which test was used

    for lag_idx = 1:nLags
        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        if length(group1) > 1 && length(group2) > 1
            % Normality check
            is_normal1 = lillietest(group1) == 0;
            is_normal2 = lillietest(group2) == 0;

            if is_normal1 && is_normal2
                % Use t-test
                [~, p, ~, stats] = ttest2(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = stats.tstat;
                test_type(lag_idx) = "ttest";
            else
                % Use non-parametric Mann–Whitney U test
                p = ranksum(group1, group2);
                p_vals(lag_idx) = p;
                t_stats(lag_idx) = NaN;  % t-stat not applicable
                test_type(lag_idx) = "ranksum";
            end
        end
    end

    % Store in structure
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;
    xcf_by_group(seed_idx).test_type = test_type;
end


sig_lags = find(p_vals <0.05);
    if ~isempty(sig_lags)
        fprintf("Seed %s: %d significant lag(s) found (p < 0.05)\n", ...
            xcf_by_group(seed_idx).seed_name, length(sig_lags));
        fprintf("  Lags: %s\n", mat2str(xcf_by_group(seed_idx).lags(sig_lags)));
        fprintf("  p-values: %s\n", mat2str(p_vals(sig_lags), 3));  % rounded to 3 decimals
    else
        fprintf("Seed %s: No significant lags (p < 0.05)\n", ...
            xcf_by_group(seed_idx).seed_name);
    end

%% run permuted p-val
% sub x lag (randomise grp IDs to determine significance)
sig_permute= nan(nSeeds,length(lags));
pval_permute= nan(nSeeds,length(lags));
for seed_idx = 1:nSeeds
    del_grp_crosscor= xcf_by_group(seed_idx).delirium;
    health_grp_crosscor= xcf_by_group(seed_idx).health;
    for i=1:length(lags)
        [sig_permute(seed_idx,i),pval_permute(seed_idx,i)]=perm_code(del_grp_crosscor(:,i),health_grp_crosscor(:,i),1000);
    end
    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;
end
%rename so variables don't get confused
xcf_diffpeakts_grp_stats = xcf_by_group;
%% Figure plot all subjects delta peak time-series
seed_names = {'LC-PPN', 'LC-nbM', 'LC-ACh', 'PPN-LC', 'nbM-LC', 'ACh-LC', 'LC+nbM'};

% plot all 7
nSubs = size(results_peakdiffts, 1);    % e.g., 119
nSeeds = size(results_peakdiffts, 2);   % e.g., 7

figure;
set(gcf, 'color', 'w');             % White background

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);        % 2x4 layout (with 1 slot left over)
    hold on;

    % Extract lags from the first subject
    lags = results_peakdiffts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(sub, :) = this_entry.mean_xcf;
            plot(lags, this_entry.mean_xcf);  
        end
    end

    % Plot group mean
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakdiffts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Optional: Adjust layout or add legend if desired
legend('Subjects', 'Group Mean', 'Location', 'best');


nSubs = size(results_peakdiffts, 1);    % e.g., 119
nSeeds = size(results_peakdiffts, 2);   % e.g., 4



figure;
set(gcf, 'color', 'w');             % White background
% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.7608, 0.2392, 0.1804;
    0.7608, 0.2392, 0.1804;
    0.403921568627451, 0.6392156862745098, 0.6235294117647059;      % Dark green for seed 2 (PPN)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
    0.11764705882352941, 0.5098039215686274, 0.49411764705882355;% Light green for seed 4 (ACh)
    0.403921568627451, 0.49411764705882355,0.6470588235294118 %pale blue
];
for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject
    lags = results_peakdiffts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(sub, :) = this_entry.mean_xcf;
            plot(lags, this_entry.mean_xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakdiffts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Add a single legend to the last subplot (optional)
subplot(2, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');

filename = fullfile(pwd, 'crosscorr_diffpeakts_allsub_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');

% simplicity avg.
figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts(1, seed_idx).lags;
    all_mean_xcf = nan(nSubs, length(lags));  % Preallocate

    for sub = 1:nSubs
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(sub, :) = this_entry.mean_xcf;
        end
    end

    % Plot group mean XCF
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    plot(lags, group_mean, 'Color', seed_colors(seed_idx, :), 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakdiffts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean XCF');
end

filename = fullfile(pwd, 'crosscorr_diffpeakts_mean_xcf.svg');  % Saves in current folder
print(gcf, filename, '-dsvg');



% plot based upon grouping
group_subs = health_sub;
nSubs_group = length(group_subs);  % Number of subjects in group
nSeeds = size(results_peakdiffts, 2);

figure;
set(gcf, 'color', 'w');  % White background

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;  % Red for seed 1 (LC)
    0.7608, 0.2392, 0.1804;
    0.7608, 0.2392, 0.1804;
    0.4039, 0.6392, 0.6235;  % Dark green for seed 2 (PPN)
    0.2627, 0.4196, 0.3725;  % Medium green for seed 3 (nbM)
    0.1176, 0.5098, 0.4941;  % Light green for seed 4 (ACh)
    0.4039, 0.4941, 0.6471;  % Pale blue
];

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject in group (to get lag range)
    lags = results_peakdiffts(group_subs(1), seed_idx).lags;
    all_mean_xcf = nan(nSubs_group, length(lags));  % Preallocate

    for i = 1:nSubs_group
        sub = group_subs(i);  % Actual subject index
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_mean_xcf(i, :) = this_entry.mean_xcf;
            plot(lags, this_entry.mean_xcf, 'Color', seed_colors(seed_idx,:), 'LineWidth', 0.5);  
        end
    end

    % Plot group mean in black on top
    group_mean = mean(all_mean_xcf, 1, 'omitnan');
    h2 = plot(lags, group_mean, 'k', 'LineWidth', 2);

    % Titles and labels
    title(strrep(results_peakdiffts(group_subs(1), seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Cross-corr');
end

% Add a single legend to the last subplot (optional)
subplot(2, 4, nSeeds);
legend(h2, 'Group Mean', 'Location', 'best');

%plotting the group means overlayed
figure;
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(2, 4, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakdiffts(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_health(i, :) = this_entry.mean_xcf;
        end
    end

    % Compute group means
    group_mean_del = mean(all_xcf_del, 1, 'omitnan');
    group_sem_del  = std(all_xcf_del, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_del),1));
    group_mean_health = mean(all_xcf_health, 1, 'omitnan');
    group_sem_health  = std(all_xcf_health, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_xcf_health),1));

        % Plot shaded error region for health
    fill([lags, fliplr(lags)], ...
         [group_mean_health + group_sem_health, fliplr(group_mean_health - group_sem_health)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot shaded error region for delirium
    fill([lags, fliplr(lags)], ...
         [group_mean_del + group_sem_del, fliplr(group_mean_del - group_sem_del)], ...
         seed_colors(seed_idx,:), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot both lines with same color but different line styles
    h1 = plot(lags, group_mean_del, '--', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % dashed
    h2 = plot(lags, group_mean_health, '-', 'Color', seed_colors(seed_idx,:), 'LineWidth', 2);  % solid

    % Titles and labels
    title(strrep(results_peakdiffts(1, seed_idx).seed_name, '_', '\_'));
    xlabel('Lag');
    ylabel('Mean XCF');

    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

% Save figure as .svg
filename = fullfile(pwd, 'crosscorr_diffpeakts_mean_xcf_del_vs_health_with_sem.svg');
print(gcf, filename, '-dsvg');





nSeeds = size(results_peakdiffts, 2);

total_peaks_matrix = zeros(nSubs, nSeeds);

for sub = 1:nSubs
    for seed = 1:nSeeds
        total_peaks_matrix(sub, seed) = results_peakdiffts(sub, seed).total_peaks;
    end
end

total_peaks_del = sum(total_peaks_matrix(delirium_sub, :), 2);
total_peaks_nondel = sum(total_peaks_matrix(health_sub, :), 2);

mean_peaks_del     = mean(total_peaks_matrix(delirium_sub, :), 1);
mean_peaks_nondel  = mean(total_peaks_matrix(health_sub, :), 1);
nondel_color = [253, 205, 154] / 255;  % #fdcd9a
del_color    = [174, 216, 230] / 255;  % #aed8e6

% Compute group stats
mean_peaks_nondel = mean(total_peaks_matrix(health_sub, :), 1);
mean_peaks_del    = mean(total_peaks_matrix(delirium_sub, :), 1);

std_peaks_nondel  = std(total_peaks_matrix(health_sub, :), 0, 1);
std_peaks_del     = std(total_peaks_matrix(delirium_sub, :), 0, 1);

group_means = [mean_peaks_nondel; mean_peaks_del];  % [2 x nSeeds]
group_stds  = [std_peaks_nondel;  std_peaks_del];   % [2 x nSeeds]

% Create bar plot
figure;
set(gcf,'color','w')
hold on;
b = bar(group_means', 'grouped');
b(1).FaceColor = nondel_color;
b(2).FaceColor = del_color;

% Calculate error bar positions
num_groups = size(group_means, 1);  % 2 groups
num_seeds  = size(group_means, 2);
group_width = min(0.8, num_groups / (num_groups + 1.5));
x = nan(num_groups, num_seeds);

for i = 1:num_groups
    x(i, :) = (1:num_seeds) - group_width/2 + (2*i-1) * group_width / (2*num_groups);
end

% Add error bars
for i = 1:num_groups
    errorbar(x(i, :), group_means(i, :), group_stds(i, :), ...
        'k.', 'LineWidth', 1, 'CapSize', 8);
end

% Aesthetics
set(gca, 'XTick', 1:num_seeds, 'XTickLabel', seed_names, 'XTickLabelRotation', 45);
ylabel('Average Number of Peaks');
legend({'Non-Delirious', 'Delirious'}, 'Location', 'northeast');
title('Average Peaks per Seed by Group');
grid on;
hold off;

p_vals = zeros(1, length(seed_names));
t_stats = zeros(1, length(seed_names));

for s = 1:length(seed_names)
    grp1 = total_peaks_matrix(health_sub, s);
    grp2 = total_peaks_matrix(delirium_sub, s);
    
    [~, p, ~, stats] = ttest2(grp1, grp2);
    
    p_vals(s) = p;
    t_stats(s) = stats.tstat;
end

% Display results
T = table(seed_names(:), t_stats(:), p_vals(:), 'VariableNames', {'Seed', 'T_stat', 'P_value'})
% Mann-Whitney Test - non-parametric assumption
p_vals_ranksum = zeros(1, length(seed_names));

for s = 1:length(seed_names)
    grp1 = total_peaks_matrix(health_sub, s);
    grp2 = total_peaks_matrix(delirium_sub, s);
    
    p = ranksum(grp1, grp2);  % Mann-Whitney U test
    p_vals_ranksum(s) = p;
end

% Display results
T_nonparam = table(seed_names(:), p_vals_ranksum(:), 'VariableNames', {'Seed', 'P_value_ranksum'})




%% Plot lall of the lags against the correlation with the BOLD activity
         % retain plots between iterations
    ts = squeeze(ts_all(sub, :, :));            % [time x region]
    ts_compare = ts(:,7:end-5);                % trim edges

    % Define seed time series (adjusting for trimmed indices)
    lc_ts  = mean(ts_compare([485,486],:), 1);
    ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
seed_ts = nbm_ts-lc_ts;
seed_ts = lc_ts-nbm_ts;

figure;              % create a new figure
hold on;    
for vv = 1:400
    [xcf, lags] = crosscorr(part(vv,:), seed_ts);

    plot(lags, xcf); % plot each curve
    [~, I] = max(abs(xcf));
    allLags(1,vv) = I;
    allCorrs(1,vv) = xcf(I);
end

hold off;
xlabel('Lags');
ylabel('Cross-correlation');
title('Cross-correlation of pc and lc - nbm ts')

save('cross_corr_analysis.mat', '-v7.3');


% After calculation of mean Participation Loop through each voxel and
% caclulating max to recreate travelling wave
allLags =nan(1,400);
for vv = 1:400
[xcf,lags] = crosscorr(mean(part_ts1),cortSig(:,vv));

%plot fig 1C
plot(xcf,lags) 
[M,I] = max(xcf);
allLags(vv) = lags(I);
end







% mean PC following peak timeseries(LC - nbM):






%% Look at the cross-correlation between the dFC of the LC, PPN and nbM (to cortex) vs PC