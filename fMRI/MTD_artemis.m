%% Artemis MTD Calculation
%% 1. Set-up Variables
addpath(genpath('project/IPOD_B3/code/'));
filepath = '/scratch/IPOD_B3/V1_MRI/timeseries/';
addpath(genpath(filepath))
addpath(genpath('/scratch/IPOD_B3/V1_MRI/timeseries/MTD/'))

%subject list - from dir
cd /scratch/IPOD_B3/V1_MRI/timeseries/
filename=dir('*timeseries.mat');

for ii=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(ii).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum)
    ses = split(2); %session ses-..
    ses = cell2mat(ses)
    
    load([subject_file]); %load in time-series; ts variable
    ts = ts(6:end -5,:); %remove first/last 5 time-points for noise
    %% MTD Coupling
    %Params for MTD
    nROI = size(ts,2);
    window=10;
    trim = 1;
    direction= 0;
    %remove first 5 and last 5 time points,remember transpose matrix to time x nodes
%      if nargin==2
%         direction=0; trim=0;
%     elseif nargin==1
%         trim=0;
%     end

    [t,nodes] = size(ts);

    %calculate temporal derivative
    td = diff(ts);

    %standardize data
    data_std = std(td);

    for i = 1:nodes
         td(:,i) = td(:,i) / data_std(1,i);
    end

    % [...] = zscore(X,FLAG,DIM) standardizes X by working along the dimension
    %    DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1 to use N.

    %functional coupling score
    fc = bsxfun(@times,permute(td,[1,3,2]),permute(td,[1,2,3]));

    %temporal smoothing (credit: T. C. O'Haver, 2008.)
    mtd_temp = zeros(nodes,nodes,t-1);
    
    for j = 1:nodes
        for k = 1:nodes
            mtd_temp(j,k,:) = smooth(squeeze(fc(:,j,k)),window);
        end
    end

    %window type (0 = middle; 1 = forward facing)
    mtd = zeros(nodes,nodes,t);

    if direction == 1
        mtd(:,:,1:t-round(window/2+1)) = mtd_temp(:,:,round(window/2+1):end);
    elseif direction == 0
        mtd(:,:,1:t-1) = mtd_temp;
    end

    %trim ends (0 = no; 1 = yes)?
    if trim == 1 && direction == 0
        mtd(:,:,t-round(window/2):end) = [];
        mtd(:,:,1:round(window/2)) = [];
    elseif trim == 1 && direction == 1
        mtd(:,:,(t-window):end) = [];
    end
    size(mtd)
    
    for nn = 1:nROI
    template = find(tril(ones(nROI))-eye(nROI));%finding unique combination of pairs
    sprintf('%d%s',nn,'completed mtd');
    end



    %flattens mtd across timepoints
    nTime = size(mtd,3);
    for tt = 1:nTime  
        temp = mtd(:,:,tt);
        mtd_flat(:,tt) = temp(template); %flattens mtd across timepoints
        sprintf('%d%s',tt,'completed mtd flat');
    end
   
    
    %save outputs
    cd /scratch/IPOD_B3/V1_MRI/timeseries/MTD/
    save([subnum '_' ses '_mtd.mat'],'mtd','-v7.3'); %save raw mtd file based on subnum
    save([subnum '_' ses '_mtd_flat.mat'],'mtd_flat','-v7.3'); %save mtd flatten across time

    clear nTime nROI nPairs temp mtd_flat mtd template
    clear nn tt
    cd .. %back to original timeseries folder
end
