%% Compute functional connectivity from CSF-regressed timeseries
% Looks for files named "*_csf4thvent_regress_timeseries.mat" in output_dir

output_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';

files = dir(fullfile(output_dir, '*_csf4thvent_regress_timeseries.mat'));
if isempty(files)
    error('No timeseries files found in %s', output_dir);
end

% Determine ROI count from first file
tmp = load(fullfile(output_dir, files(1).name));
if isfield(tmp, 'cleaned_timeseries')
    ts0 = tmp.cleaned_timeseries;
elseif isfield(tmp, 'timeseries')
    ts0 = tmp.timeseries;
else
    error('No timeseries variable found in %s', files(1).name);
end

n_rois = size(ts0, 2);
fc_all_csf4thvent_regress = zeros(length(files), n_rois, n_rois);
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

    % Functional connectivity (Pearson correlation)
    ts_corr = corr(ts);
    fc_all_csf4thvent_regress(ii, :, :) = ts_corr;
end

save(fullfile(output_dir, 'fc_all_csf4thvent_regress.mat'), 'fc_all_csf4thvent_regress', 'sub_ids');
fprintf('Saved FC for %d subjects to %s\n', length(files), fullfile(output_dir, 'fc_all_csf4thvent_regress.mat'));


%-------------------------------------------------------------------------%
%% Acompcorr Timeseries
%-------------------------------------------------------------------------%
output_dir = '/Users/ntaylor/Desktop/AASDelirium_timeseries_csf4thvent';

files = dir(fullfile(output_dir, '*_acompcor_timeseries.mat'));
if isempty(files)
    error('No timeseries files found in %s', output_dir);
end

% Determine ROI count from first file
tmp = load(fullfile(output_dir, files(1).name));
if isfield(tmp, 'ts')
    ts0 = tmp.ts;
elseif isfield(tmp, 'cleaned_timeseries')
    ts0 = tmp.cleaned_timeseries;
else
    error('No timeseries variable found in %s', files(1).name);
end

n_rois = size(ts0, 2);
fc_all_acompcor = zeros(length(files), n_rois, n_rois);
sub_ids = strings(length(files), 1);

for ii = 1:length(files)
    subject_file = files(ii).name;
    split = strsplit(subject_file, '_');
    sub_ids(ii) = string(split{1});

    data = load(fullfile(output_dir, subject_file));
    if isfield(data, 'cleaned_timeseries')
        ts = data.cleaned_timeseries;
    elseif isfield(data, 'ts')
        ts = data.ts;
    else
        error('No timeseries variable found in %s', subject_file);
    end

    % Remove first/last 5 timepoints for noise
    if size(ts, 1) > 10
        ts = ts(6:end-5, :);
    end

    % Functional connectivity (Pearson correlation)
    ts_corr = corr(ts);
    fc_all_acompcor(ii, :, :) = ts_corr;
end

output_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';
save(fullfile(output_dir, 'fc_all_acompcor.mat'), 'fc_all_acompcor', 'sub_ids');
fprintf('Saved FC for %d subjects to %s\n', length(files), fullfile(output_dir, 'fc_all_acompcor.mat'));
