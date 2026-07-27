%% Prepare Data for running 4th Ventricle CSF Extraction

%data directory
cd('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\timeseries\schaef_400')
filename=dir('*timeseries.mat');

% Initialize sub_ids as a cell array to store character arrays (subject IDs)
sub_ids = cell(length(filename), 1); 

for ii = 1:length(filename)
    % Get 'sub-...' from the name
    subject_file = filename(ii).name;  % filename of time-series data
    
    % Split the filename by underscores
    split = strsplit(subject_file, '_');
    
    % Extract 'sub-...' and convert to character array
    subnum = split{1};  % 'sub-...' section
    sub_ids{ii} = subnum;  % Store in sub_ids cell array
    
    % Extract 'ses-...' (session) and convert to character array
    ses = split{2};  % 'ses-...' section
    ses = cell2mat(ses);  % (Optional) Convert to character array if needed
end




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