%Output:     
%              ci       time-resolved community assignment
%              q        time-resolved modularity
%              p        time-resolved participation coefficient
%              z        time-resolved module-degree z-score
%              hc       cartographic profile
%              f        flexibility

cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\MTD'

%% load data MTD
filename=dir('*mtd.mat');
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum);
    ses = split(2); %session ses-..
    ses = cell2mat(ses);
    
    load([subject_file]); %load in mtd

%data size is 502 x 502 x 210
    gamma = 1;
    beta = 0.75;

    [ci,q,part,modz,hc,f] = integration_plus5(mtd,gamma,beta);

    cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\Graph_Theory\schaef_400\'
    save([subnum '_' ses '_part.mat'],"part");
    save([subnum '_' ses '_flex.mat'],"f");
    save([subnum '_' ses '_modz.mat'],"modz");
    save([subnum '_' ses '_modularity.mat'],"q");
    save([subnum '_' ses '_community.mat'],"ci");
    save([subnum '_' ses '_cartographic.mat'],"hc");

    cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\Graph_Theory\schaef_400\Additional_V1_2'
end
% additional V1 Analysis for MTD
cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\MTD\Additional_V1\'
filename=dir('*mtd.mat');
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum);
    ses = split(2); %session ses-..
    ses = cell2mat(ses);
    
    load([subject_file]); %load in mtd

%data size is 502 x 502 x 210
    gamma = 1;
    beta = 0.75;

    [ci,q,part,modz,hc,f] = integration_plus5(mtd,gamma,beta);

    cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\Graph_Theory'
    save([subnum '_' ses '_part.mat'],"part");
    save([subnum '_' ses '_flex.mat'],"f");
    save([subnum '_' ses '_modz.mat'],"modz");
    save([subnum '_' ses '_modularity.mat'],"q");
    save([subnum '_' ses '_community.mat'],"ci");
    save([subnum '_' ses '_cartographic.mat'],"hc");

    cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\MTD\Additional_V1\'
end

%% flexibility measure
filename=dir('*flex.mat');
filename([59 93])=[];
flex_all = zeros(length(filename),502);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data

    load([subject_file]); %load in mtd

flex_all(i,:)= f;
end
%to directory
cd 'C:\Users\natas\Documents\PhD\Rob_Sanders\Data\Analysis\Graph_Theory'
%% participation coefficient
cd schaef_400
filename=dir('*part.mat');
%remove shorter subject windows - sub-107 & sub-164
filename([59 93])=[];
part_all = zeros(length(filename),502,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

part_all(i,:,:)= part;
end
%% modularity measure
filename=dir('*modularity.mat');
filename([59 93])=[];
modularity_all = zeros(length(filename),210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

    modularity_all(i,:)= q;
end

%% module z measure
filename=dir('*modz.mat');
filename([59 93])=[];
modz_all = zeros(length(filename),502,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

    modz_all(i,:,:)= modz;
end

%% community measure
filename=dir('*community.mat');
filename([59 93])=[];
community_all = zeros(length(filename),502,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

    community_all(i,:,:)= ci;
end

%% cartographic measure
filename=dir('*cartographic.mat');
filename([59 93])=[];
cartographic_all = zeros(length(filename),101,101,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd

    cartographic_all(i,:,:,:)= hc;
end

