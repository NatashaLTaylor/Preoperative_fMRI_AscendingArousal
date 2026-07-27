%% Setup paths
timeseries_dir = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400';
csf_dir        = '/Users/ntaylor/Desktop/AASDelirium_timeseries_csf4thvent';
output_dir     = '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/timeseries/schaef_400/timeseries_csf4thvent_regress_brainstem';

%% Load subject IDs
% Assumes sub_ids is a cell array like {'sub-001','sub-002',...}
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/subject_data.mat');  
sub_ids=subjectV1MRIdata.sub_id;
%% Specify regions to regress CSF from
% Example: indices of regions (columns in timeseries)
%489 L/490 R - PPN
% 491 L / 492 R - nbM
% 485 L / 486 R - LC
regions_to_clean = [491:492, 485:486];  % <-- just nbM and LC regresses

%% Loop through subjects
for i = 1:length(sub_ids)
    
    sub = sub_ids(i);
    
    fprintf('Processing %s...\n', sub);
    
    %% Load ROI timeseries
    ts_file = fullfile(timeseries_dir, ...
        sprintf('%s_ses-1_timeseries.mat', sub));
    
    ts_data = load(ts_file);
    
    % Assume variable inside is called 'timeseries'
    timeseries = ts_data.timeseries;  % [T x N]
    
    %% Load CSF 4th ventricle timeseries
    csf_file = fullfile(csf_dir, ...
        sprintf('%s_4thvent_timeseries.txt', sub));
    
    csf = load(csf_file);  % [T x 1]
    
    %% Ensure matching length
    if size(timeseries,1) ~= length(csf)
        error('Timepoints do not match for %s', sub);
    end
    
    %% Add intercept to regression
    X = [csf, ones(length(csf),1)];  % [T x 2]
    
    %% Copy original timeseries
    cleaned_timeseries = timeseries;
    
    %% Regress CSF from selected regions
    for r = regions_to_clean
        
        y = timeseries(:, r);  % region signal
        
        % Solve linear regression: y = X * beta
        beta = X \ y;
        
        % Remove CSF contribution (but keep intercept)
        y_clean = y - X(:,1) * beta(1);
        
        cleaned_timeseries(:, r) = y_clean;
    end
    
    %% Save output
    output_file = fullfile(output_dir, ...
        sprintf('%s_csf4thvent_regress_timeseries.mat', sub));
    
    save(output_file, 'cleaned_timeseries');
    
end

fprintf('Done!\n');