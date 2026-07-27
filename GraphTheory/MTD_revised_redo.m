%% MTD Calculation - Revised Regression timeseries

 
% === Settings ===
timeseries_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';
output_dir = "/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/csf_regress_redo";
window = 10; % <-- set your window length
trim = 1;
direction= 0;
files = dir(fullfile(timeseries_dir, '*_csf4thvent_regress_timeseries.mat'));

% === Load first file to size arrays ===
tmp = load(fullfile(timeseries_dir, files(1).name));
if isfield(tmp, 'cleaned_timeseries')
    ts0 = tmp.cleaned_timeseries;
elseif isfield(tmp, 'timeseries')
    ts0 = tmp.timeseries;
else
    error('No timeseries variable found in %s', files(1).name);
end

% Trim
if size(ts0, 1) <= 10
    error('Not enough timepoints in first file after trimming.');
end
ts0 = ts0(6:end-5, :);

nROI = size(ts0, 2);
time = 210;

% Template for flattening (lower triangle, no diagonal)
template = find(tril(ones(nROI)) - eye(nROI));
nEdges = numel(template);

% Preallocate
nSub = length(files);
mtd_flat_all = zeros(nSub, nEdges, time, 'single'); % single to save space
sub_ids = strings(nSub, 1);

for ii = 1:nSub
    time_series_name = fullfile(timeseries_dir, files(ii).name);

    % Subject ID from filename
    split = strsplit(files(ii).name, '_');
    subnum = split{1};
    sub_ids(ii) = string(subnum);

    % Load data
    data = load(time_series_name);
    if isfield(data, 'cleaned_timeseries')
        ts = data.cleaned_timeseries;
    elseif isfield(data, 'timeseries')
        ts = data.timeseries;
    else
        error('No timeseries variable found in %s', files(ii).name);
    end

    % Remove first and last 5 timepoints
    if size(ts, 1) <= 10
        warning('Skipping %s: not enough timepoints after trimming.', files(ii).name);
        continue;
    end
    ts = ts(6:end-5, :);

    % Dimensions
    nROI = size(ts, 2);

    % Compute MTD (time x nodes assumed in coupling function)
    mtd = coupling(ts, window,direction,trim);
    time= size(mtd,3);
    % Unique lower triangle indices (exclude diagonal)
    template = find(tril(ones(nROI)) - eye(nROI));

    % Save raw MTD
    %save(fullfile(output_dir, [subnum '_mtd.mat']), 'mtd');

    % Flatten MTD across timepoints
    mtd_flat = zeros(numel(template), time);
    for tt = 1:time
        temp = mtd(:, :, tt);
        mtd_flat(:, tt) = temp(template);
        sprintf('%d%s', tt, ' completed mtd flat'); %#ok<SPRINTFN>
    end
    %save into one file
    mtd_flat_all(ii,:,:)= mtd_flat;
    % Save flattened MTD
    save(fullfile(output_dir, [subnum 'csf4thvent_mtd_flat.mat']), 'mtd_flat');

    fprintf('Finished %s\n', subnum);
end

save(fullfile(output_dir, 'mtd_flat_all_subjects_csf4thvent.mat'), 'mtd_flat_all', 'sub_ids', 'template', '-v7.3');