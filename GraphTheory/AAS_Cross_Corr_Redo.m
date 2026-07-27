% redo - 3/07/2026

% Load in relevant data
load('/Users/ntaylor/OneDrive - The University of Sydney (Staff)/Postdoc_Rob/Analysis/Demographic_Data/V1_MRI_all_demos.mat');
sub_remove_clinical = [42,43,59,63,74,93,103,104];
subjectV1MRIdata_graphtheory_orig=subjectV1MRIdata;
subjectV1MRIdata_graphtheory_orig(sub_remove_clinical,:)=[];
a= find(subjectV1MRIdata_graphtheory_orig.sub_id==49);
b= find(subjectV1MRIdata_graphtheory_orig.sub_id==50);
sub_remove_grp1=[a,b];
subjectV1MRIdata_graphtheory_update = subjectV1MRIdata_graphtheory_orig;
subjectV1MRIdata_graphtheory_update(sub_remove_grp1,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata_graphtheory_update(1:size(subjectV1MRIdata_graphtheory_update,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1);
health_sub_update =bin_delirium_all~=1; 
health_sub = find(health_sub_update==1);


% Original Timeseries:
% load in the original timeseries for subject analysed -
% dir: '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/'




% Updated csf regress Timeseries:
% Assume data is stored in a matrix 'data' of size T x N
output_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';

% files = dir(fullfile(output_dir, '*_csf4thvent_regress_timeseries.mat'));
% if isempty(files)
%     error('No timeseries files found in %s', output_dir);
% end
% 
% max_timepoints = 221;
% num_rois = 502;  % number of ROIs
% ts_all = nan(length(files), num_rois, max_timepoints);  % preallocate with NaNs
% 
% sub_ids = strings(length(files), 1);
% 
% for ii = 1:length(files)
%     subject_file = files(ii).name;
%     split = strsplit(subject_file, '_');
%     sub_ids(ii) = string(split{1});
% 
%     data = load(fullfile(output_dir, subject_file));
%     if isfield(data, 'cleaned_timeseries')
%         ts = data.cleaned_timeseries;
%     elseif isfield(data, 'timeseries')
%         ts = data.timeseries;
%     else
%         error('No timeseries variable found in %s', subject_file);
%     end
% 
%     % Remove first/last 5 timepoints for noise
%     if size(ts, 1) > 10
%         ts = ts(6:end-5, :);
%     end
% 
%     % Check length
%     if size(ts, 1) == max_timepoints
%         ts_all(ii, :, :) = ts';  % transpose to roi x time
%     else
%         ts_all(ii, :, :) = NaN;  % entire time series replaced with NaNs
%     end
% end
% 
%sub_remove_ts = [sub_remove_grp1,88]; %remove sub-049, sub-050, sub-164
%ts_all(sub_remove_ts,:,:)=[];

% load in the updated csf timeseries:
load("/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/csf_regress_redo/ts_csf4thvent_regress_subs.mat")

load('participation_coeff_csf4thvent_all.mat')
%part_all = 117 x 502 x 210




% ============================ %
%% Peak Normal BOLD in Timeseries
% ============================ %
%1.  peaks in BOLD timeseries
nSubs=size(part_all,1);
seed_names = {'LC','nbM'};
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
std_thresh=2; % vary the threshold

%BOLD trimming ts_compare = ts(:,7:end-5);   

for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,7:end-5);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);

    seed_ts_list = {lc_ts, nbm_ts};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s}; %select time-series

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
        [~,peak_inds] =  findpeaks(seed_ts,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(seed_ts)); %find peaks of acc
        %total peaks
        sum_pks = length(peak_inds);
        xcf_accum = [];
        std_ts = std(ts);
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
        results_peakts(sub, s).subject   = sub;
        results_peakts(sub, s).seed_name = seed_names{s};
        results_peakts(sub, s).timeseries = seed_ts;
        results_peakts(sub, s).peak_inds = peak_inds;
        results_peakts(sub, s).total_peaks = sum_pks;
        results_peakts(sub, s).std_ts    = std_ts;
        results_peakts(sub, s).xcf  = xcf_accum;
        results_peakts(sub, s).mean_xcf  = mean_xcf;
        results_peakts(sub, s).lags      = lags;
        results_peakts(sub, s).max_corr  = maxCorr;
        results_peakts(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

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

    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);

    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);

    for lag_idx = 1:nLags

        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        % Only run test if enough observations
        if length(group1) > 1 && length(group2) > 1

            [~, p, ~, stats] = ttest2(group1, group2);

            p_vals(lag_idx) = p;
            t_stats(lag_idx) = stats.tstat;

        end
    end

    % Store results
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;

end

% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end

%rename so variables don't get confused
xcf_peakts_grp_stats = xcf_by_group;


% plot

%plotting the group means overlayed
figure('Theme','light');
set(gcf, 'color', 'w');
set(gcf, 'Color', 'w');
seed_names = {'LC', 'nbM'};

for s = 1:2
    subplot(1, 2, s);
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

nSubs = size(results_peakts, 1);    % e.g., 117
nSeeds = size(results_peakts, 2);   % e.g., 4

figure('Theme','light');
set(gcf, 'color', 'w');   
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
];

for seed_idx = 1:nSeeds
    subplot(1, 2, seed_idx);
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
    xticks(lags)
    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

filename = fullfile(pwd, 'crosscorr_peakts_2_stdthresh_mean_xcf_del_vs_health_wtih_sem_updated_csfregress.svg');
print(gcf, filename, '-dsvg');

%% Mediation Analysis- variables
% extract variables for mediation analysis
nSubs = size(part_all, 1);
seed_names = {'LC', 'nbM'};

% Preallocate structure array
% NEW: added peak_frequency, std_seed, mean_peak_amplitude fields
results_peakts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [], ...
    'peak_inds', [], ...
    'total_peaks', [], ...
    'peak_freq', [], ...       % NEW: peaks per unit time (normalised)
    'std_seed', [], ...             % FIX: renamed from std_ts, now seed-specific
    'mean_peak_amp', [], ...  % NEW: mean amplitude of suprathreshold peaks
    'xcf', [], ...
    'mean_xcf', [], ...
    'lags', [], ...
    'max_corr', [], ...
    'max_lag', [], ...
    'mean_pc',[]);

window_size = 5;
std_thresh = 2;

for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));
    ts_compare = ts(:, 7:end-5);
    part = squeeze(part_all(sub, :, :));
    mean_pc = mean(part, 1);

    % Seed time series
    lc_ts  = mean(ts_compare([485,486], :), 1);
    nbm_ts = mean(ts_compare([491,492], :), 1);

    seed_ts_list = {lc_ts, nbm_ts};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % FIX: standard deviation of THIS seed's time series
        seed_std = std(seed_ts);

        % Detect peaks above 2 SDs from mean
        [peak_vals, peak_inds] = findpeaks(seed_ts, ...
            'MinPeakDistance', 1, ...
            'MinPeakHeight', std_thresh * seed_std);
        % NOTE: peak_vals now captured (first output of findpeaks)

        % Total peak count
        sum_pks = length(peak_inds);

        % NEW: Peak frequency — peaks normalised by time series length
        % This accounts for any differences in usable time points across subjects
        n_timepoints = length(seed_ts);
        peak_freq = sum_pks / n_timepoints;

        % NEW: Mean amplitude of suprathreshold peaks
        % peak_vals already contains only the amplitudes at detected peaks
        if sum_pks > 0
            mean_peak_amp = mean(peak_vals);
        else
            mean_peak_amp = NaN;
        end

        % Cross-correlation around peaks (unchanged)
        xcf_accum = [];

        for i = 1:length(peak_inds)
            idx = peak_inds(i);
            start_idx = idx - window_size;
            end_idx   = idx + window_size;

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
            maxLag  = lags(I);
        else
            mean_xcf = nan(1, 2*window_size+1);
            lags = -window_size:window_size;
            maxCorr = NaN;
            maxLag  = NaN;
        end

        % Save all results
        results_peakts(sub, s).subject             = sub;
        results_peakts(sub, s).seed_name           = seed_names{s};
        results_peakts(sub, s).timeseries          = seed_ts;
        results_peakts(sub, s).peak_inds           = peak_inds;
        results_peakts(sub, s).total_peaks         = sum_pks;
        results_peakts(sub, s).peak_freq      = peak_freq;      % NEW
        results_peakts(sub, s).std_seed             = seed_std;       % FIX
        results_peakts(sub, s).mean_peak_amp = mean_peak_amp;  % NEW
        results_peakts(sub, s).xcf                 = xcf_accum;
        results_peakts(sub, s).mean_xcf            = mean_xcf;
        results_peakts(sub, s).lags                = lags;
        results_peakts(sub, s).max_corr            = maxCorr;
        results_peakts(sub, s).max_lag             = maxLag;
        results_peakts(sub, s).mean_pc       = mean(mean_pc);  
    end

    fprintf("Subject finished %d\n", sub);
end



%% NEW: Extract mediation variables into a table for export
% =========================================================

%% Export mediation variables to CSV
% Build table for export

export_rows = [];
varNames = {'subject', 'seed_name', 'delirium', 'mean_pc', ...
            'peak_freq', 'std_seed', 'mean_peak_amp', 'total_peaks', ...
            'max_corr', 'max_lag'};

for sub = 1:nSubs
    for s = 1:length(seed_names)
        r = results_peakts(sub, s);
        row = {r.subject, ...
               r.seed_name, ...
               bin_delirium_all(sub), ...
               r.mean_pc, ...
               r.peak_freq, ...
               r.std_seed, ...
               r.mean_peak_amp, ...
               r.total_peaks, ...
               r.max_corr, ...
               r.max_lag};
        export_rows = [export_rows; row];
    end
end

T = cell2table(export_rows, 'VariableNames', varNames);
writetable(T, 'mediation_variables_crosscorr.csv');
fprintf('Exported mediation_variables.csv with %d rows\n', height(T));


%% redo mediation variable extract
% just take the time-lag +1 PC values in correspondance to the peak in BOLD
results_peakts(nSubs, length(seed_names)) = struct( ...
    'subject', [], ...
    'seed_name', '', ...
    'timeseries', [], ...
    'peak_inds', [], ...
    'total_peaks', [], ...
    'peak_freq', [], ...
    'std_seed', [], ...
    'mean_peak_amp', [], ...
    'mean_pc_post_peak', [], ...   % NEW: mean PC at 1 step after each peak
    'xcf', [], ...
    'mean_xcf', [], ...
    'lags', [], ...
    'max_corr', [], ...
    'max_lag', [], ...
    'mean_pc', []);

window_size = 5;
std_thresh = 2;

for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));
    ts_compare = ts(:, 7:end-5);
    part = squeeze(part_all(sub, :, :));
    mean_pc = mean(part, 1);

    % Seed time series
    lc_ts  = mean(ts_compare([485,486], :), 1);
    nbm_ts = mean(ts_compare([491,492], :), 1);

    seed_ts_list = {lc_ts, nbm_ts};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % FIX: standard deviation of THIS seed's time series
        seed_std = std(seed_ts);

        % Detect peaks above 2 SDs from mean
        [peak_vals, peak_inds] = findpeaks(seed_ts, ...
            'MinPeakDistance', 1, ...
            'MinPeakHeight', std_thresh * seed_std);
        % NOTE: peak_vals now captured (first output of findpeaks)

        % Total peak count
        sum_pks = length(peak_inds);
        % Extract participation coefficient at 1 timepoint after each peak
        post_peak_pc = [];
        for i = 1:length(peak_inds)
            idx = peak_inds(i) + 1;  % one step ahead
            if idx <= length(mean_pc)
                post_peak_pc = [post_peak_pc, mean_pc(idx)];
            end
        end

        if ~isempty(post_peak_pc)
            mean_pc_post_peak = mean(post_peak_pc);
        else
            mean_pc_post_peak = NaN;
        end
        % NEW: Peak frequency — peaks normalised by time series length
        % This accounts for any differences in usable time points across subjects
        n_timepoints = length(seed_ts);
        peak_freq = sum_pks / n_timepoints;

        % NEW: Mean amplitude of suprathreshold peaks
        % peak_vals already contains only the amplitudes at detected peaks
        if sum_pks > 0
            mean_peak_amp = mean(peak_vals);
        else
            mean_peak_amp = NaN;
        end

        % Cross-correlation around peaks (unchanged)
        xcf_accum = [];

        for i = 1:length(peak_inds)
            idx = peak_inds(i);
            start_idx = idx - window_size;
            end_idx   = idx + window_size;

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
            maxLag  = lags(I);
        else
            mean_xcf = nan(1, 2*window_size+1);
            lags = -window_size:window_size;
            maxCorr = NaN;
            maxLag  = NaN;
        end

        % Save all results
        results_peakts(sub, s).subject             = sub;
        results_peakts(sub, s).seed_name           = seed_names{s};
        results_peakts(sub, s).timeseries          = seed_ts;
        results_peakts(sub, s).peak_inds           = peak_inds;
        results_peakts(sub, s).total_peaks         = sum_pks;
        results_peakts(sub, s).peak_freq      = peak_freq;      % NEW
        results_peakts(sub, s).std_seed             = seed_std;       % FIX
        results_peakts(sub, s).mean_peak_amp = mean_peak_amp;  % NEW
        results_peakts(sub, s).xcf                 = xcf_accum;
        results_peakts(sub, s).mean_xcf            = mean_xcf;
        results_peakts(sub, s).lags                = lags;
        results_peakts(sub, s).max_corr            = maxCorr;
        results_peakts(sub, s).max_lag             = maxLag;
        results_peakts(sub, s).mean_pc       = mean(mean_pc); 
        results_peakts(sub, s).mean_pc_post_peak = mean_pc_post_peak;  % NEW
    end

    fprintf("Subject finished %d\n", sub);
end


%% Export mediation variables to CSV
varNames = {'subject', 'seed_name', 'delirium', 'mean_pc', 'mean_pc_post_peak', ...
            'peak_freq', 'std_seed', 'mean_peak_amp', 'total_peaks', ...
            'max_corr', 'max_lag'};

export_rows = [];

for sub = 1:nSubs
    for s = 1:length(seed_names)
        r = results_peakts(sub, s);
        row = {r.subject, ...
               r.seed_name, ...
               bin_delirium_all(sub), ...
               r.mean_pc, ...
               r.mean_pc_post_peak, ...    % NEW
               r.peak_freq, ...
               r.std_seed, ...
               r.mean_peak_amp, ...
               r.total_peaks, ...
               r.max_corr, ...
               r.max_lag};
        export_rows = [export_rows; row];
    end
end

T = cell2table(export_rows, 'VariableNames', varNames);
writetable(T, 'mediation_variables_postpeak_pc.csv');
fprintf('Exported mediation_variables.csv with %d rows\n', height(T));




% replaces for the PC variables:

T2 = T;
% replace pc with the average of this significantly different ones:

% average of the significantly different PCs:

% If data is 502 x 210
repeated_sig_pc = repmat(sig_pc_time, 1, 1, 117);  % gives 502 x 210 x 117

% Then permute to get 117 x 502 x 210
repeated_sig_pc  = permute(repeated_sig_pc , [3 1 2]);
part_all_sig_time = part_all.*repeated_sig_pc;

part_avg_all_sig = mean(part_all_sig_time,3);
avg_part_avg_all_sig = mean(part_avg_all_sig,2);
T2.mean_pc = repelem(avg_part_avg_all_sig, 2);
writetable(T, 'mediation_variables_crosscorr_sig_pc_all.csv');

% take the average of just timepoints that are not 0
part_all_sig_time(part_all_sig_time ==0)=NaN;
part_avg_all_sig = mean(part_all_sig_time,3,'omitnan');
avg_part_avg_all_sig = mean(part_avg_all_sig,2,'omitnan');
% remove the 0-lag correlation by calculating a partial cross-correlation:

T2.mean_pc = repelem(avg_part_avg_all_sig, 2);
writetable(T, 'mediation_variables_crosscorr_sig_pc.csv');


%2.  peaks in delta BOLD timeseries
seed_names = {'LC-nbM', 'nbM-LC', 'LC+nbM'};

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
std_thresh=2; % standard deviation threshold
%calculating across each sub
for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,6:end-6);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    %ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    ach_ts = mean(ts_compare([489:492],:), 1);

    %lcppn_ts = lc_ts - ppn_ts;
    lcnbm_ts = lc_ts - nbm_ts;
    lcach_ts = lc_ts - ach_ts;
    %ppnlc_ts = ppn_ts - lc_ts;
    nbmlc_ts = nbm_ts - lc_ts;
    achlc_ts = ach_ts - lc_ts;
    lcnbm_combine = lc_ts + nbm_ts;

    %seed_ts_list = {lcppn_ts, lcnbm_ts, lcach_ts, ppnlc_ts, nbmlc_ts, achlc_ts, lcnbm_combine};

    seed_ts_list = {lcnbm_ts, nbmlc_ts,lcnbm_combine};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
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


% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end


% figure plotting
group_subs = health_sub;
nSubs_group = length(group_subs);  % Number of subjects in group
nSeeds = size(results_peakdiffts, 2);

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;  % Red for seed 1 (LC)
    0.2627, 0.4196, 0.3725;  % Medium green for seed 3 (nbM)
    0.4039, 0.4941, 0.6471;  % Pale blue
];

%plotting the group means overlayed
figure('Theme','light');
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, nSeeds, seed_idx);
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
    xticks(lags);
    if seed_idx == nSeeds
        legend([h2 h1], {'Healthy', 'Delirium'}, 'Location', 'best');
    end
end

% Save figure as .svg
filename = fullfile(pwd, 'crosscorr_diffpeakts_2sd_mean_xcf_del_vs_health_with_sem_updated.svg');
print(gcf, filename, '-dsvg');





% ============================ %
%% Peak Second Derivative BOLD in Timeseries
% ============================ %

% peaks in the second derivative of the BOLD signals
%1.  peaks in BOLD timeseries
nSubs=size(part_all,1);
seed_names = {'LC','nbM'};
% Preallocate structure array
results_peakts_deriv(nSubs, length(seed_names)) = struct( ...
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
std_thresh=2.5; % vary the threshold

for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,7:end-5);
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    
    % second derivative:
    %sigAcc = gradient(gradient(sig));
    lc_ts_sigacc = gradient(gradient(lc_ts));
    nbm_ts_sigacc = gradient(gradient(nbm_ts));

    seed_ts_list = {lc_ts_sigacc, nbm_ts_sigacc};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s}; %select time-series

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
        [~,peak_inds] =  findpeaks(seed_ts,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(seed_ts)); %find peaks of acc
        %total peaks
        sum_pks = length(peak_inds);
        xcf_accum = [];
        std_ts = std(ts);
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
        results_peakts_deriv(sub, s).subject   = sub;
        results_peakts_deriv(sub, s).seed_name = seed_names{s};
        results_peakts_deriv(sub, s).timeseries = seed_ts;
        results_peakts_deriv(sub, s).peak_inds = peak_inds;
        results_peakts_deriv(sub, s).total_peaks = sum_pks;
        results_peakts_deriv(sub, s).std_ts    = std_ts;
        results_peakts_deriv(sub, s).xcf  = xcf_accum;
        results_peakts_deriv(sub, s).mean_xcf  = mean_xcf;
        results_peakts_deriv(sub, s).lags      = lags;
        results_peakts_deriv(sub, s).max_corr  = maxCorr;
        results_peakts_deriv(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts_deriv(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
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

    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);

    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);

    for lag_idx = 1:nLags

        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        % Only run test if enough observations
        if length(group1) > 1 && length(group2) > 1

            [~, p, ~, stats] = ttest2(group1, group2);

            p_vals(lag_idx) = p;
            t_stats(lag_idx) = stats.tstat;

        end
    end

    % Store results
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;

end

% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end

%rename so variables don't get confused
xcf_peakts_grp_stats = xcf_by_group;

% figure generations

figure('Theme','light');
seed_names = {'LC', 'nbM'};
set(gcf, 'color', 'w');             % White background
% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
];


% plot

figure;
set(gcf, 'Color', 'w');
seed_names = {'LC', 'nbM'};

for s = 1:2
    subplot(1, 2, s);
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

nSubs = size(results_peakts, 1);    % e.g., 119
nSeeds = size(results_peakts, 2);   % e.g., 4

figure('Theme','light');
set(gcf, 'color', 'w');   
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
];

for seed_idx = 1:nSeeds
    subplot(1, 2, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts_deriv(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
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

filename = fullfile(pwd, 'crosscorr_peakts_2ndderiv_1_5stdthresh_mean_xcf_del_vs_health_wtih_sem_updated_csfregress.svg');
print(gcf, filename, '-dsvg');




%2.  peaks in delta BOLD timeseries
seed_names = {'LC-nbM', 'nbM-LC', 'LC+nbM'};

% Preallocate structure array
results_peakdiffts_deriv(nSubs, length(seed_names)) = struct( ...
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
std_thresh=2; % standard deviation threshold
%calculating across each sub
for sub = 1:nSubs
    ts = squeeze(ts_all(sub, :, :));         
    ts_compare = ts(:,6:end-6);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    %ach_ts = mean(ts_compare([489:492],:), 1);

    lcnbm_ts = lc_ts - nbm_ts;
    nbmlc_ts = nbm_ts - lc_ts;
    lcnbm_combine = lc_ts + nbm_ts;
    
    % second derivative:
    %sigAcc = gradient(gradient(sig));
    lc_ts_sigacc = gradient(gradient(lc_ts));
    nbm_ts_sigacc = gradient(gradient(nbm_ts));
    lcnbmcombine_ts_sigacc = gradient(gradient(lcnbm_combine));

    seed_ts_list = {lc_ts_sigacc, nbm_ts_sigacc, lcnbmcombine_ts_sigacc};


    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
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
        results_peakdiffts_deriv(sub, s).subject   = sub;
        results_peakdiffts_deriv(sub, s).seed_name = seed_names{s};
        results_peakdiffts_deriv(sub, s).timeseries = seed_ts;
        results_peakdiffts_deriv(sub, s).peak_inds = peak_inds;
        results_peakdiffts_deriv(sub, s).total_peaks = sum_pks;
        results_peakdiffts_deriv(sub, s).std_ts    = std_ts;
        results_peakdiffts_deriv(sub, s).xcf  = xcf_accum;
        results_peakdiffts_deriv(sub, s).mean_xcf  = mean_xcf;
        results_peakdiffts_deriv(sub, s).lags      = lags;
        results_peakdiffts_deriv(sub, s).max_corr  = maxCorr;
        results_peakdiffts_deriv(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

% is there a significant difference between the groups

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts_deriv(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakdiffts_deriv(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakdiffts_deriv(sub, seed_idx);
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


% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end


% figure plotting
group_subs = health_sub;
nSubs_group = length(group_subs);  % Number of subjects in group
nSeeds = size(results_peakdiffts, 2);

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;  % Red for seed 1 (LC)
    0.2627, 0.4196, 0.3725;  % Medium green for seed 3 (nbM)
    0.4039, 0.4941, 0.6471;  % Pale blue
];

%plotting the group means overlayed
figure('Theme','light');
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, nSeeds, seed_idx);
    hold on;

    % Extract lags from the first subject (assumed same across all)
    lags = results_peakdiffts_deriv(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakdiffts_deriv(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakdiffts_deriv(sub, seed_idx);
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
filename = fullfile(pwd, 'crosscorr_diffpeakts_2ndderiv_2sd_mean_xcf_del_vs_health_with_sem_updated.svg');
print(gcf, filename, '-dsvg');







%% Original Timeseries Loaded:
 % original timeseries for nbM - as there aren't any major predictions of
 % CSF 4thvent affecting it's timeseries but still use the recalculated
 % participation coefficient

 cd '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400'
filename = dir("*timeseries.mat");

max_timepoints = 221;
num_rois = 502;  % number of ROIs
ts_all_orig = nan(length(filename), num_rois, max_timepoints);  % preallocate with NaNs
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
        ts_all_orig(ii, :, :) = ts';  % transpose to roi x time
    else
        ts_all_orig(ii, :, :) = NaN;  % entire time series replaced with NaNs
    end
end
%now remove relevant IDs
sub_remove_clinical = [30,31,42,43,59,63,74,93,103,104];
% and 30/31 needs to be removed
ts_all_orig(sub_remove_clinical,:,:)=[];

% ============================ %
%% Peak Normal BOLD in Timeseries
% ============================ %
%1.  peaks in BOLD timeseries
nSubs=size(part_all,1);
seed_names = {'LC','nbM'};
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
std_thresh=2.5; % vary the threshold

for sub = 1:nSubs
    ts = squeeze(ts_all_orig(sub, :, :));         
    ts_compare = ts(:,7:end-5);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);

    seed_ts_list = {lc_ts, nbm_ts};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s}; %select time-series

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
        [~,peak_inds] =  findpeaks(seed_ts,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(seed_ts)); %find peaks of acc
        %total peaks
        sum_pks = length(peak_inds);
        xcf_accum = [];
        std_ts = std(ts);
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
        results_peakts(sub, s).subject   = sub;
        results_peakts(sub, s).seed_name = seed_names{s};
        results_peakts(sub, s).timeseries = seed_ts;
        results_peakts(sub, s).peak_inds = peak_inds;
        results_peakts(sub, s).total_peaks = sum_pks;
        results_peakts(sub, s).std_ts    = std_ts;
        results_peakts(sub, s).xcf  = xcf_accum;
        results_peakts(sub, s).mean_xcf  = mean_xcf;
        results_peakts(sub, s).lags      = lags;
        results_peakts(sub, s).max_corr  = maxCorr;
        results_peakts(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

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

    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);

    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);

    for lag_idx = 1:nLags

        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        % Only run test if enough observations
        if length(group1) > 1 && length(group2) > 1

            [~, p, ~, stats] = ttest2(group1, group2);

            p_vals(lag_idx) = p;
            t_stats(lag_idx) = stats.tstat;

        end
    end

    % Store results
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;

end

% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end

%rename so variables don't get confused
xcf_peakts_grp_stats = xcf_by_group;


% plot

%plotting the group means overlayed
seed_names = {'LC', 'nbM'};

figure('Theme','light');
set(gcf, 'color', 'w');   
seed_colors = [
    0.7608, 0.2392, 0.1804;      % Red for seed 1 (LC)
    0.2627450980392157, 0.4196078431372549, 0.37254901960784315;      % Medium green for seed 3 (nbM)
];

for seed_idx = 1:nSeeds
    subplot(1, 2, seed_idx);
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

filename = fullfile(pwd, 'crosscorr_peakts_2_5stdthresh_mean_xcf_del_vs_health_wtih_sem_updated_csfregress.svg');
print(gcf, filename, '-dsvg');








%2.  peaks in delta BOLD timeseries
seed_names = {'LC-nbM', 'nbM-LC', 'LC+nbM'};

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
std_thresh=1; % standard deviation threshold
%calculating across each sub
for sub = 1:nSubs
    ts = squeeze(ts_all_orig(sub, :, :));         
    ts_compare = ts(:,6:end-6);              
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    %ppn_ts = mean(ts_compare([489,490],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    ach_ts = mean(ts_compare([489:492],:), 1);

    %lcppn_ts = lc_ts - ppn_ts;
    lcnbm_ts = lc_ts - nbm_ts;
    lcach_ts = lc_ts - ach_ts;
    %ppnlc_ts = ppn_ts - lc_ts;
    nbmlc_ts = nbm_ts - lc_ts;
    achlc_ts = ach_ts - lc_ts;
    lcnbm_combine = lc_ts + nbm_ts;

    %seed_ts_list = {lcppn_ts, lcnbm_ts, lcach_ts, ppnlc_ts, nbmlc_ts, achlc_ts, lcnbm_combine};

    seed_ts_list = {lcnbm_ts, nbmlc_ts,lcnbm_combine};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s};

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
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


% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end


% figure plotting
group_subs = health_sub;
nSubs_group = length(group_subs);  % Number of subjects in group
nSeeds = size(results_peakdiffts, 2);

% Define custom colors for seeds
seed_colors = [
    0.7608, 0.2392, 0.1804;  % Red for seed 1 (LC)
    0.2627, 0.4196, 0.3725;  % Medium green for seed 3 (nbM)
    0.4039, 0.4941, 0.6471;  % Pale blue
];

%plotting the group means overlayed
figure('Theme','light');
set(gcf, 'color', 'w');  % White background

for seed_idx = 1:nSeeds
    subplot(1, nSeeds, seed_idx);
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
filename = fullfile(pwd, 'crosscorr_diffpeakts_1p5sd_mean_xcf_del_vs_health_with_sem_updated.svg');
print(gcf, filename, '-dsvg');





% ============================ %
%% Peak Second Derivative BOLD in Timeseries
% ============================ %

% peaks in the second derivative of the BOLD signals
%1.  peaks in BOLD timeseries
nSubs=size(part_all,1);
seed_names = {'LC','nbM'};
% Preallocate structure array
results_peakts_deriv(nSubs, length(seed_names)) = struct( ...
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
std_thresh=1; % vary the threshold

for sub = 1:nSubs
    ts = squeeze(ts_all_orig(sub, :, :));         
    ts_compare = ts(:,7:end-5);
    part = squeeze(part_all(sub,:,:));       
    mean_pc = mean(part, 1);                 

    % Seed time series
    lc_ts  = mean(ts_compare([485,486],:), 1);
    nbm_ts = mean(ts_compare([491,492],:), 1);
    
    % second derivative:
    %sigAcc = gradient(gradient(sig));
    lc_ts_sigacc = gradient(gradient(lc_ts));
    nbm_ts_sigacc = gradient(gradient(nbm_ts));

    seed_ts_list = {lc_ts_sigacc, nbm_ts_sigacc};

    for s = 1:length(seed_ts_list)
        seed_ts = seed_ts_list{s}; %select time-series

        % Detect peaks above 2 SDs from mean
        stdThresh =std_thresh;
        [~,peak_inds] =  findpeaks(seed_ts,'MinPeakDistance',1,'MinPeakheight',stdThresh*std(seed_ts)); %find peaks of acc
        %total peaks
        sum_pks = length(peak_inds);
        xcf_accum = [];
        std_ts = std(ts);
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
        results_peakts_deriv(sub, s).subject   = sub;
        results_peakts_deriv(sub, s).seed_name = seed_names{s};
        results_peakts_deriv(sub, s).timeseries = seed_ts;
        results_peakts_deriv(sub, s).peak_inds = peak_inds;
        results_peakts_deriv(sub, s).total_peaks = sum_pks;
        results_peakts_deriv(sub, s).std_ts    = std_ts;
        results_peakts_deriv(sub, s).xcf  = xcf_accum;
        results_peakts_deriv(sub, s).mean_xcf  = mean_xcf;
        results_peakts_deriv(sub, s).lags      = lags;
        results_peakts_deriv(sub, s).max_corr  = maxCorr;
        results_peakts_deriv(sub, s).max_lag   = maxLag;
    end

    fprintf("Subject finished %d\n", sub);
end

xcf_by_group = struct();

nSeeds = length(seed_names);
for seed_idx = 1:nSeeds
    % Extract lags from the first subject (assumed same across all)
    lags = results_peakts_deriv(1, seed_idx).lags;

    % Preallocate
    all_xcf_del = nan(length(delirium_sub), length(lags));
    all_xcf_health = nan(length(health_sub), length(lags));

    % Delirium subjects
    for i = 1:length(delirium_sub)
        sub = delirium_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
        if ~isempty(this_entry.mean_xcf)
            all_xcf_del(i, :) = this_entry.mean_xcf;
        end
    end

    % Healthy subjects
    for i = 1:length(health_sub)
        sub = health_sub(i);
        this_entry = results_peakts_deriv(sub, seed_idx);
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

    del = xcf_by_group(seed_idx).delirium;
    health = xcf_by_group(seed_idx).health;

    nLags = length(xcf_by_group(seed_idx).lags);

    p_vals = nan(1, nLags);
    t_stats = nan(1, nLags);

    for lag_idx = 1:nLags

        group1 = del(:, lag_idx);
        group2 = health(:, lag_idx);

        % Remove NaNs
        group1 = group1(~isnan(group1));
        group2 = group2(~isnan(group2));

        % Only run test if enough observations
        if length(group1) > 1 && length(group2) > 1

            [~, p, ~, stats] = ttest2(group1, group2);

            p_vals(lag_idx) = p;
            t_stats(lag_idx) = stats.tstat;

        end
    end

    % Store results
    xcf_by_group(seed_idx).p_vals = p_vals;
    xcf_by_group(seed_idx).t_stats = t_stats;

end

% prints out the significance 
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


% print out sig for permutation
alpha=0.05;
for seed_idx = 1:nSeeds
    
    del_grp_crosscor = xcf_by_group(seed_idx).delirium;
    health_grp_crosscor = xcf_by_group(seed_idx).health;

    for i = 1:length(lags)
        [sig_permute(seed_idx,i), pval_permute(seed_idx,i)] = ...
            perm_code(del_grp_crosscor(:,i), health_grp_crosscor(:,i), 1000);
    end

    xcf_by_group(seed_idx).permuted_sig = sig_permute;
    xcf_by_group(seed_idx).permuted_pval = pval_permute;

    % ---- PRINT SIGNIFICANT LAGS ----
    sig_lag_idx = find(pval_permute(seed_idx,:) < alpha);
    sig_lags = lags(sig_lag_idx);

    fprintf('\nSeed %d significant lags:\n', seed_idx);
    disp(sig_lags);
end

%rename so variables don't get confused
xcf_peakts_grp_stats = xcf_by_group;
