%% Graph Theory Analysis

%% dFC Analysis

cd '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400'

%load in clinical data
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Demographic_Data/V1_MRI_all_demos.mat');

% just look at within DMN nodes
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
sub_remove(3,:)=[]; % need to remove 59, as that has already been removed
subjectV1MRIdata(sub_remove,:)=[];

sub_remove_clinical = [42,43,59,63,74,93,103,104];
filename = dir('*mtd_flat.mat');
%filename([59,93])=[]; %remove as the time length not the same
mtd_all_flat =zeros(size(filename,1),125751,210);
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    load([subject_file]); %load in mtd 502 x 502 x 210
    mtd_all_flat(i,:,:,:)= mtd_flat;
    sprintf('%d',i)
end
%% dFC analysis 
avg_dFC = mean(mtd_all_flat,3); %avg across time
avg_dFC([sub_remove,:])=[];
avg_dFC_del = mean(mtd_all_flat(delirium_sub,:,:),3);
avg_avg_dFC_del = matify((mean(avg_dFC_del)'),502);
avg_dFC_health = mean(mtd_all_flat(health_sub,:,:),3);
avg_avg_dFC_health = matify((mean(avg_dFC_health)'),502);
std_dFC = std(mtd_all_flat,[],3);
std_dFC([sub_remove],:)=[];
std_dFC_del = std_dFC(delirium_sub,:);
avg_std_dFC_del = matify(mean(std_dFC_del)',502);
std_dFC_health = std_dFC(health_sub,:); 
avg_std_dFC_health = matify(mean(std_dFC_health)',502);
%stats
for i=1:size(avg_dFC,2)
   [sig_dFC_grp(i,:),pval_dFC_grp(i,:)]=perm_code(avg_dFC_del(:,i),avg_dFC_health(:,i),1000);
end

for i=1:size(std_dFC,2)
   [sig_std_dFC_grp(i,:),pval_std_dFC_grp(i,:)]=perm_code(std_dFC_del(:,i),std_dFC_health(:,i),1000);
end

a = mean(avg_dFC_del);
sig_avg_dFC_del = matify((a'.*sig_dFC_grp),502);
b = mean(avg_dFC_health);
sig_avg_dFC_health = matify((b'.*sig_dFC_grp),502);
diff_sig_avg_dFC = sig_avg_dFC_del-sig_avg_dFC_health;

a = mean(std_dFC_del);
sig_std_dFC_del = matify((a'.*sig_std_dFC_grp),502);
b = mean(std_dFC_health);
sig_std_dFC_health = matify((b'.*sig_std_dFC_grp),502);
diff_sig_std_dFC = sig_std_dFC_del-sig_std_dFC_health;

%2. Model for dFC
%demographs. data to compare
clinical_tbl = table(subjectV1MRIdata.patient_age,subjectV1MRIdata.drs_total_calculated,subjectV1MRIdata.overall_peak_drs_calc2,subjectV1MRIdata.new_nsqip_d, ...
    'VariableNames',{'age','drstotal_baseline','peak_drs_postop','nsqip_d'});


dFC_avg_tbl = array2table(avg_dFC); %avg. across time
mdl_coeff_avgdFC = zeros(size(avg_dFC,2),4);
for i=1:size(avg_dFC,2)
    %load in individual FC edges for each sub
    roi_tbl = [clinical_tbl dFC_avg_tbl(:,i)];
    mdl= fitglm(roi_tbl,'ResponseVar',"peak_drs_postop",'Distribution','normal');
    mdl_coeff_avgdFC(i,1) = mdl.Coefficients.Estimate(5); %beta coeffs
    mdl_coeff_avgdFC(i,2) = mdl.Coefficients.pValue(5); %pvalue
    mdl_coeff_avgdFC(i,3) = mdl.Coefficients.SE(5); %SE
    mdl_coeff_avgdFC(i,4) = mdl.Coefficients.tStat(5); %tStat
end

%sig dFC from model
sig_mdl_avgdFC = double(mdl_coeff_avgdFC(:,2)<0.05);
sig_mdl_coef_avgdFC = sig_mdl_avgdFC.*mdl_coeff_avgdFC(:,1);
sig_mdl_coef_avgdFC_mat =matify(sig_mdl_coef_avgdFC,502);
[h, crit_p,~,adj_p]=fdr_bh(mdl_coeff_avgdFC(:,2),0.05,'pdep','yes');  %FDR correction for multiple comparisons


%sig. SD from model
dFC_std_tbl = array2table(std_dFC); %sd across time
mdl_coeff_stddFC = zeros(size(std_dFC,2),4);
for i=1:size(std_dFC,2)
    %load in individual FC edges for each sub
    roi_tbl = [clinical_tbl dFC_std_tbl(:,i)];
    mdl= fitglm(roi_tbl,'ResponseVar',"peak_drs_postop",'Distribution','normal');
    mdl_coeff_stddFC(i,1) = mdl.Coefficients.Estimate(5); %beta coeffs
    mdl_coeff_stddFC(i,2) = mdl.Coefficients.pValue(5); %pvalue
    mdl_coeff_stddFC(i,3) = mdl.Coefficients.SE(5); %SE
    mdl_coeff_stddFC(i,4) = mdl.Coefficients.tStat(5); %tStat
end
sig_mdl_stddFC = double(mdl_coeff_stddFC(:,2)<0.05);
sig_mdl_coef_stddFC = sig_mdl_stddFC.*mdl_coeff_stddFC(:,1);
sig_mdl_coef_stddFC_mat =matify(sig_mdl_coef_stddFC,502);
[h, crit_p,~,adj_p]=fdr_bh(mdl_coeff_stddFC(:,2),0.05,'pdep','yes');  %FDR correction for multiple comparisons

%% Ascending Arousal System
for i=1:119
    sig_std_dFC_all(i,:) = sig_std_dFC_grp'.*std_dFC(i,:);
end

for i=1:119
    a(i,:,:) = matify(sig_std_dFC_all(i,:)',502);
end

for i=1:119
    a(i,:,:)=matify(std_dFC(i,:,:)',502);
end

sig_aas_std_dFC = a(:,:,483:end);
flat_aas_std_dFC = reshape(sig_aas_std_dFC,119,502*20);

[rho_aas_std_dFC,pval_aas_std_dFC]=corr(flat_aas_std_dFC,peak_drs_overall,'rows','complete');

sig_rho_aas_std_dFC = rho_aas_std_dFC.*(double(pval_aas_std_dFC<0.05));
reshape_sig_rho_aas_std_dFC= reshape(sig_rho_aas_std_dFC,502,20);
%figure
figure
set(gcf,'color','w')
imagesc(reshape_sig_rho_aas_std_dFC(voltron_order(:,3),:))
ylabel('ROIs')
title('Sig. Rho SD dFC AAS vs peak drs')
colormap(Color_Salmon_Blue)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)



% focus on LC/nbM & PPN dFC
% Define seed region index (change as needed)
%489 L/490 R - PPN
% 491 L / 492 R - nbM
% 485 L / 486 R - LC



for i=1:size(avg_dFC,2)
   [sig_dFC_grp(i,:),pval_dFC_grp(i,:)]=perm_code(avg_dFC_del(:,i),avg_dFC_health(:,i),1000);
end

for i=1:size(std_dFC,2)
   [sig_std_dFC_grp(i,:),pval_std_dFC_grp(i,:)]=perm_code(std_dFC_del(:,i),std_dFC_health(:,i),1000);
end




%% focus on DMN foucs for dFC

%1. Look across all L/R nodes of DMN to rest of the brain
DMN_ids = double(find(voltron_order(:,2)==10)); %make sure getting correct IDs for left and right, did it wrong on Friday 24/05
%Left Default Network A/C: 149-194 or
%281*502 = 120982
a = zeros(119,502,502);
for i=1:119
    a(i,:,:)=matify(std_dFC(i,:)',502);
end
mat_std_dFC = a;
DMN_dFC_std = mat_std_dFC(:,:,[149:194,358:390]); %just take DMN 77 nodes

DMN_dfc_std_del = DMN_dFC_std(delirium_sub,:,:); %this contains the 0s for within self
DMN_dfc_std_health = DMN_dFC_std(health_sub,:,:);
diff_DMN_dfc_std = DMN_dfc_avg_del - DMN_dfc_avg_health;

%stats
%flatten
flat_DMN_dfc_std_del = reshape(DMN_dfc_std_del,31,502*size(DMN_dfc_std_health,3));
flat_DMN_dfc_std_health = reshape(DMN_dfc_std_del,31,502*size(DMN_dfc_std_health,3));

for i=1:size(DMN_dfc_std_del,2)
   for k=1:size(DMN_dfc_std_del,3)
   [sig_DMN_dFC_std_grp(i,k),pval_DMN_dFC_std_grp(i,k)]=perm_code(DMN_dfc_std_del(:,i,k),DMN_dfc_std_health(:,i,k),1000);
   end
end

for i=1:size(flat_DMN_dfc_std_del,2)
       [sig_DMN_dFC_std_grp_2(i,:),pval_DMN_dFC_std_grp_2(i,:)]=perm_code(flat_DMN_dfc_std_del(:,i),flat_DMN_dfc_std_health(:,i),1000);
end

%the significant nodes from the original permutation across all 502x502
sig_diff_DMN_dfc_std_1 =sig_std_dFC_del(:,[149:194,358:390])- sig_std_dFC_health(:,[149:194,358:390]);

%correlation with peak drs calc.
load('clinical_data.mat')
flat_DMN_dFC_std = reshape(DMN_dFC_std,119,502*79);

[rho_dfc_std_DMN,pval_dfc_std_DMN]=corr(flat_DMN_dFC_std,peak_drs_overall,'rows','complete');

sig_rho_dfc_std_DMN = rho_dfc_std_DMN.*(double(pval_dfc_std_DMN<0.05));
reshape_sig_rho_dfc_std_DMN = reshape(sig_rho_dfc_std_DMN,502,79);

% avg. dFC first for DMN nodes altogether
avg_DMN_dFC_std = mean(DMN_dFC_std,3);
avg_DMN_dFC_std_del = avg_DMN_dFC_std(delirium_sub,:);
avg_DMN_dFC_std_health = avg_DMN_dFC_std(health_sub,:);

for i=1:size(avg_DMN_dFC_std_del,2)
    [sig_avg_DMN_dFC_std(i,:),pval_avg_DMN_dFC_std(i,:)]=perm_code(avg_DMN_dFC_std_del(:,i),avg_DMN_dFC_std_health(:,i),1000);
end

figure
set(gcf,'color','w')
subplot(1,2,1)
a = mean(avg_DMN_dFC_std_del,1)'.*sig_avg_DMN_dFC_std;
imagesc(a(voltron_order(:,3),:))
ylabel('ROIs')
title('Sig. Avg. DMN SD dFC Delirium')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
subplot(1,2,2)
b = mean(avg_DMN_dFC_std_health,1)'.*sig_avg_DMN_dFC_std;
imagesc(a(voltron_order(:,3),:))
ylabel('ROIs')
title('Sig. Avg. DMN SD dFC Health')
colormap(Color_Blues)


subplot(1,3,3)
diff = a-b;
imagesc(diff(voltron_order(:,3),:))
ylabel('ROIs')
title('Sig. Diff Avg. DMN SD dFC Del - Health')
colormap(Color_Salmon_Blue)



figure
set(gcf,'color','w')
subplot(1,3,1)
imagesc(DMN_dfc_avg_del(voltron_order(:,3),:))
xlabel('All DMN nodes')
ylabel('ROIs')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum



%figure
figure
set(gcf,'color','w')
subplot(1,2,1)
%imagesc(avg_std_dFC_del(voltron_order(:,3),voltron_order(:,3)));
%imagesc(sig_avg_dFC_del(voltron_order(:,3),voltron_order(:,3)));
imagesc(sig_std_dFC_del(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('SD. dFC delirium')
%title('Avg. dFC delirium')
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
%imagesc(avg_std_dFC_health(voltron_order(:,3),voltron_order(:,3)));
%imagesc(sig_avg_dFC_health(voltron_order(:,3),voltron_order(:,3)));
imagesc(sig_std_dFC_health(voltron_order(:,3),voltron_order(:,3)));
xlabel('ROIs')
ylabel('ROIs')
title('SD. dFC healthy')
%title('Avg. dFC healthy')
colormap(CustomColormap4)
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

% look at dFC of ascending arousal system structures
figure
set(gcf,'color','w')
subplot(1,3,1)
%imagesc(sig_avg_dFC_del(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
imagesc(sig_std_dFC_del(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
ylabel('All ROIs')
xlabel('AAS ROIs')
title('Sig. SD. dFC Delirium')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([2.5,2.5],[1,502],'Color','black','LineWidth',1) %DR
line([4.5,4.5],[1,502],'Color','black','LineWidth',1) %LC
line([6.5,6.5],[1,502],'Color','black','LineWidth',1) %VTA
line([8.5,8.5],[1,502],'Color','black','LineWidth',1) %PPN
line([10.5,10.5],[1,502],'Color','black','LineWidth',1) %nbM
line([12.5,12.5],[1,502],'Color','black','LineWidth',1) %PO
line([14.5,14.5],[1,502],'Color','black','LineWidth',1) %PAG
line([16.5,16.5],[1,502],'Color','black','LineWidth',1) % Hypothalamus
line([18.5,18.5],[1,502],'Color','black','LineWidth',1) % STN
%last is rednucleus
subplot(1,3,2)
imagesc(sig_std_dFC_health(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
%imagesc(sig_avg_dFC_health(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
ylabel('All ROIs')
xlabel('AAS ROIs')
title('Sig. SD. dFC Health')
%title('Sig. Avg. dFC Health')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([2.5,2.5],[1,502],'Color','black','LineWidth',1) %DR
line([4.5,4.5],[1,502],'Color','black','LineWidth',1) %LC
line([6.5,6.5],[1,502],'Color','black','LineWidth',1) %VTA
line([8.5,8.5],[1,502],'Color','black','LineWidth',1) %PPN
line([10.5,10.5],[1,502],'Color','black','LineWidth',1) %nbM
line([12.5,12.5],[1,502],'Color','black','LineWidth',1) %PO
line([14.5,14.5],[1,502],'Color','black','LineWidth',1) %PAG
line([16.5,16.5],[1,502],'Color','black','LineWidth',1) % Hypothalamus
line([18.5,18.5],[1,502],'Color','black','LineWidth',1) % STN
subplot(1,3,3)
imagesc(diff_sig_std_dFC(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
%imagesc(diff_sig_avg_dFC(voltron_order(:,3),voltron_order(483:end,3))) %just ascending arousal
ylabel('All ROIs')
xlabel('AAS ROIs')
title('Sig. Diff SD, dFC Del - Health')
%title('Sig. Diff Avg. dFC Del - Health')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([2.5,2.5],[1,502],'Color','black','LineWidth',1) %DR
line([4.5,4.5],[1,502],'Color','black','LineWidth',1) %LC
line([6.5,6.5],[1,502],'Color','black','LineWidth',1) %VTA
line([8.5,8.5],[1,502],'Color','black','LineWidth',1) %PPN
line([10.5,10.5],[1,502],'Color','black','LineWidth',1) %nbM
line([12.5,12.5],[1,502],'Color','black','LineWidth',1) %PO
line([14.5,14.5],[1,502],'Color','black','LineWidth',1) %PAG
line([16.5,16.5],[1,502],'Color','black','LineWidth',1) % Hypothalamus
line([18.5,18.5],[1,502],'Color','black','LineWidth',1) % STN
colormap(CustomColormap4)






%% clinical data
% load in clinical data
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Demographic_Data\V1_MRI_all_demos.mat');
%need to remove subs 107,164
subjectV1MRIdata([59,93],:)=[];
%find IDs that have 'NaN'
sub_remove = find(subjectV1MRIdata.overall_delirious_ever=='NA');
subjectV1MRIdata(sub_remove,:)=[];
bin_delirium_all=table2array(subjectV1MRIdata(1:size(subjectV1MRIdata,1),"bin_delirium"));
delirium_sub = find(bin_delirium_all==1); %1 if delirious ever
health_sub =bin_delirium_all~=1; 
health_sub = find(health_sub==1);
%sub-groups of delirium
hyperactive_ids =[2;5;13;15;25;26;34;38;49;62;92;138];
hyperactive_sub = find(ismember(subjectV1MRIdata.sub_id,hyperactive_ids)==1); %location in file for subjects
hypoactive_sub = find(~ismember(delirium_sub,hyperactive_sub)==1); %remainder assume hypoactive

peak_drs_overall = subjectV1MRIdata.overall_peak_drs_calc2;

%list of subject IDs to check clinical data
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum);
    split = strsplit(subnum,'-');
    sub_all{i}=cell2mat(split(2));
end
sub_all = sub_all';



%% Participation Coefficient Analysis
cd('C:\Users\natas\OneDrive - The University of Sydney (Staff)\Postdoc_Rob\Analysis\Graph_Theory\schaef_400\')
load('participation_coeff_all.mat') %loads in part_all
part_all(sub_remove,:,:)=[]; %remove subjects without clinical data
%group separation
part_hyper = part_all(hyperactive_sub,:,:);
part_hypo = part_all(hypoactive_sub,:,:);
part_del = part_all(delirium_sub,:,:);
part_healthy = part_all(health_sub,:,:);

max =max(mean(part_all(delirium_sub,:,:)),[],'all')
min =min(mean(part_all(hypoactive_sub,:,:)),[],'all')
clim = [min max]


% 1. Group Difference - just delirium vs. healthy
%avg. PC between grps
load("colormaps.mat")
figure
set(gcf,'Color','w')
subplot(1,2,1)
imagesc(squeeze(mean(part_del).*hh_pc_del));
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC delirium subs')
subplot(1,2,2)
imagesc(squeeze(mean(part_healthy).*hh_pc_del));
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC healthy sub')
colormap(CustomColormap3)
% plot of permutation significance
%load in network connectivity order
load('C:\Users\natas\OneDrive - The University of Sydney (Staff)\PhD\Code\Parcellations\voltron_ordered.mat')
figure
set(gcf,'Color','w')
subplot(1,2,1)
a= squeeze(mean(part_del));
imagesc(a(voltron_order(:,3),:).*sig_pc_time(voltron_order(:,3),:));
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC delirium subs')
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
%line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
%line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
%line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
%line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
%line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
%line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
subplot(1,2,2)
b= squeeze(mean(part_healthy));
imagesc(b(voltron_order(:,3),:).*sig_pc_time(voltron_order(:,3),:));
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC healthy sub')
colormap(CustomColormap3)

%2. variability of PC between grps

std_pc_del = std(part_del,[],3); %variability of PC across time for each ROI
std_pc_health = std(part_healthy,[],3);
load('colormaps2.mat')
figure
set(gcf,'Color','w')
subplot(1,2,1)
imagesc(squeeze(mean(std_pc_del))');
xlabel('Time')
ylabel('Participation Coefficient')
title('SD PC delirium subs')
colormap(Color_Salmon_Blue2)
subplot(1,2,2)
imagesc(squeeze(mean(std_pc_health))');
xlabel('Time')
ylabel('Participation Coefficient')
title('SD PC healthy sub')
colormap(Color_Salmon_Blue2)

load('schaef_order.mat')

figure
set(gcf,'Color','w')
subplot(1,3,1)
a= squeeze(mean(std_pc_del))'.*sig_pc_sd;
imagesc(a(voltron_order(:,3)));
xlabel('Time')
ylabel('Participation Coefficient')
title('SD PC delirium subs')
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)
subplot(1,3,2)
b=squeeze(mean(std_pc_health))'.*sig_pc_sd;
imagesc(b(voltron_order(:,3)));
xlabel('Time')
ylabel('Participation Coefficient')
title('SD PC healthy sub')
colormap(Color_Salmon_Blue2)
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)
subplot(1,3,3)
imagesc(diff_sd_pc(voltron_order(:,3)));
xlabel('Time')
ylabel('Participation Coefficient')
title('SD PC del - healthy sub')
hold on
line([0,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([0,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([0,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([0,502],[244,244],'Color','black','LineWidth',1) % limbic
line([0,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([0,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([0,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([0,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,502],[430,430],'Color','black','LineWidth',1) %thal
line([0,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([0,502],[482,482],'Color','black','LineWidth',1)




a = squeeze(mean(std_pc_del)').*sig_pc_sd;
RB_surf_schaef(a(1:400),'sig_variability_pc_delirium')
colormap(CustomColormap4)

b = squeeze(mean(std_pc_health)').*sig_pc_sd;
RB_surf_schaef(b(1:400),'sig_variability_pc_healthy')
colormap(CustomColormap4)

diff_sd_pc = a - b; %del - healthy


%stats 
[hh_pc_del,pp_pc_del]=ttest2(part_del,part_healthy);
[hh_std_pc_del,pp_std_pc_del]=ttest2(std_pc_del,std_pc_health);

% group-wise permutation significant of PC edges
nROIs = size(part_all,2);
for i=1:nROIs
   for k=1:size(part_all,3)
   [sig_pc_time(i,k),pval_pc_time(i,k)]=perm_code(part_del(:,i,k),part_healthy(:,i,k),1000);
   end
end
%permute sd of pc
for i=1:nROIs
   [sig_pc_sd(i,:),pval_pc_sd(i,:)]=perm_code(std_pc_del(:,i),std_pc_health(:,i),1000);
end

% 3.  pearson correlation
reshape_all_part = reshape(part_all,107,502*210);
%avg. pc per ROI
reshape_avg_part = mean(part_all,3);
[rho_avg_pc_drs_total,pval_avg_pc_drs_total]=corr(reshape_avg_part,drs_total_calculated2,'rows','complete');
[rho_pc_drs_total,pval_pc_drs_total]=corr(reshape_all_part,drs_total_calculated2,'rows','complete');


[rho_avg_pc_drs_total,pval_avg_pc_drs_total]=corr(reshape_avg_part,drs_total_calculated2,'rows','complete');
figure
set(gcf,'Color','w');
imagesc(reshape(rho_pc_drs_total,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS')
colormap(CustomColormap)

%% Linear model analysis
avg_part_all = mean(part_all,3);
std_part_all = std(part_all,[],3);
a = avg_part_all;
a = std_part_all;
columns_to_remove = all(a == 1);
a(:, columns_to_remove) = []; %120 x10020
b = 1:size(a,2);
c = vertcat(b,a); %need additional column read it into it
writematrix(c, 'avg_pc_subjects.csv');
writematrix(c, 'std_pc_subjects.csv');
writetable(subjectV1MRIdata,'subject_data_dFC.csv'); %dFC specific data into table

%load file from LM in R
std_pc_lm_model = readmatrix("permute_lm_std_pc_coeffs_pval.csv");

sig_permute = double(std_pc_lm_model(:,2)<0.05);
sig_coefs = sig_permute.*std_pc_lm_model(:,1);
sig_coefs(length(std_pc_lm_model),:)=[];

figure
set(gcf,'color','w')
imagesc(sig_coefs(voltron_order(:,3),:))
title('Sig. Permuted Beta SD PC vs delirium severity')
hold on
%line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,1.5],[47,47],'Color','black','LineWidth',1)
line([0,1.5],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([0,1.5],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([0,1.5],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([0,1.5],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([0,1.5],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([0,1.5],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([0,1.5],[400,400],'Color','black','LineWidth',1) %temp par net
%line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([0,1.5],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
%line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,1.5],[430,430],'Color','black','LineWidth',1) %thal
%line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([0,1.5],[454,454],'Color','black','LineWidth',1) %basal gang.
%line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([0,1.5],[482,482],'Color','black','LineWidth',1) %cerebellum
%line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(Color_Salmon_Blue)
% brain plot
RB_surf_schaef(sig_coefs(1:400,:),'')
subcort_plot(sig_coefs)

%% Module Degree Z
load('modz_all.mat')
modz_all(sub_remove,:,:)=[];
avg_modz_all = mean(modz_all,3); %avg. across time
std_modz_all = std(modz_all,[],3); %variability across time
a = avg_modz_all;
a = std_modz_all;
b = 1:size(a,2);
c = vertcat(b,a); %need additional column read it into i
writematrix(c, 'std_modz_subjects.csv');

%load avg. module degree lm beta coefs
avg_modz_lm_model = readmatrix("permute_lm_avg_modz_coeffs_pval.csv");

sig_permute = double(avg_modz_lm_model(:,2)<0.05);
sig_coefs = sig_permute.*avg_modz_lm_model(:,1);
sig_coefs(length(avg_modz_lm_model),:)=[];

%load std (variability) module degree lm beta coefs
std_modz_lm_model = readmatrix("permute_lm_std_modz_coeffs_pval.csv");

sig_permute = double(std_modz_lm_model(:,2)<0.05);
sig_coefs_std_modz = sig_permute.*std_modz_lm_model(:,1);
sig_coefs_std_modz(length(std_modz_lm_model),:)=[];



figure
set(gcf,'color','w')
imagesc(sig_coefs_std_modz(voltron_order(:,3),:))
title('Sig. Permuted Beta SD ModZ vs delirium severity')
hold on
%line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,1.5],[47,47],'Color','black','LineWidth',1)
line([0,1.5],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([0,1.5],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([0,1.5],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([0,1.5],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([0,1.5],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([0,1.5],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([0,1.5],[400,400],'Color','black','LineWidth',1) %temp par net
%line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([0,1.5],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
%line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,1.5],[430,430],'Color','black','LineWidth',1) %thal
%line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([0,1.5],[454,454],'Color','black','LineWidth',1) %basal gang.
%line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([0,1.5],[482,482],'Color','black','LineWidth',1) %cerebellum
%line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(Color_Salmon_Blue)
xticklabels([])
% brain plot
RB_surf_schaef(sig_coefs_std_modz(1:400,:),'')
subcort_plot(sig_coefs_std_modz)

%% Flexibility Calculation
%load in the file
load('flexibility_all.mat')
flex_all(sub_remove,:)=[]; %becomes 119 x 502
b = 1:size(flex_all,2);
c = vertcat(b,flex_all); %need additional column read it into i
writematrix(c, 'flex_all_subjects.csv');
%load in the params. from the permuted lm
%load lm beta coefs
flex_lm_model = readmatrix("permute_lm_flex_coeffs_pval.csv");

sig_permute = double(flex_lm_model(:,2)<0.05);
locs = sig_permute=1;
sig_coefs_flex = sig_permute.*flex_lm_model(:,1);
sig_coefs_flex(length(flex_lm_model),:)=[];

figure
set(gcf,'color','w')
imagesc(sig_coefs_flex(voltron_order(:,3),:))
title('Sig. Permuted Beta Flexibility vs delirium severity')
hold on
%line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([0,1.5],[47,47],'Color','black','LineWidth',1)
line([0,1.5],[117,117],'Color','black','LineWidth',1) %somat mot a/b
%line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([0,1.5],[169,169],'Color','black','LineWidth',1) %dorsatten
%line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([0,1.5],[220,220],'Color','black','LineWidth',1) %sal ventatten
%line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([0,1.5],[244,244],'Color','black','LineWidth',1) % limbic
%line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([0,1.5],[305,305],'Color','black','LineWidth',1) %all Cont A-C
%line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([0,1.5],[384,384],'Color','black','LineWidth',1) %all default nets
%line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
line([0,1.5],[400,400],'Color','black','LineWidth',1) %temp par net
%line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([0,1.5],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
%line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([0,1.5],[430,430],'Color','black','LineWidth',1) %thal
%line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([0,1.5],[454,454],'Color','black','LineWidth',1) %basal gang.
%line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([0,1.5],[482,482],'Color','black','LineWidth',1) %cerebellum
%line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(Color_Salmon_Blue)
xticklabels([])
% brain plot
RB_surf_schaef(sig_coefs_flex(1:400,:),'')
subcort_plot(sig_coefs_flex)






%% Figures
% pval only significant
sig_rho_pc_drs_total = rho_pc_drs_total.*(pval_pc_drs_total<0.05);
figure
set(gcf,'Color','w');
imagesc(reshape(sig_rho_pc_drs_total,502,210));
xlabel('Time')
ylabel('ROIs')
title('Sig. Correlation between PC to total DRS')
colormap(CustomColormap)




% 2. Group difference (sub-types of delirium 

%group separation of delirium types
figure
set(gcf,'Color','w');
subplot(1,4,1)
imagesc(squeeze(mean(part_all)),clim);
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC all subs')
subplot(1,4,2)
imagesc(squeeze(mean(part_all(hyperactive_sub,:,:))),clim);
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC hyperactive subs')
subplot(1,4,3)
imagesc(squeeze(mean(part_all(hypoactive_sub,:,:))),clim);
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC hypoactive subs')
subplot(1,4,4)
imagesc(squeeze(mean(part_all(health_sub,:,:))),clim);
xlabel('Time')
ylabel('Participation Coefficient')
title('Avg. PC healthy subs')

%stats

[hh_pc_del,pp_pc_del]=ttest2(part_hypo,part_healthy);


% group-wise statistical difference

for i=1:502
   for k=1:210
   [sig_pc_hypo(i,k),pval_pc_hypo(i,k)]=perm_code(part_hypo(:,i,k),part_healthy(:,i,k),1000);
   end
end

max = max(squeeze(mean(part_hypo)).*sig_pc_hypo,[],'all');
a = sig_pc_hypo.*squeeze(mean(part_hypo)-mean(part_healthy));
min = min(a,[],'all');

figure
clim = [0 max];
clim1 = [min max];
set(gcf,'color','w')
s(1)=subplot(1,3,1)
imagesc(sig_pc_hypo.*squeeze(mean(part_hypo)),clim)
title('hypoactive delirium')
colormap(s(1),CustomPastelBlue)
xlabel('TRs')
ylabel('ROIs')
s(2)=subplot(1,3,2)
imagesc(sig_pc_hypo.*squeeze(mean(part_healthy)),clim)
title('healthy')
colormap(s(2),CustomPastelBlue)
xlabel('TRs')
ylabel('ROIs')
colorbar
s(3) =subplot(1,3,3)
imagesc(sig_pc_hypo.*squeeze(mean(part_hypo)-mean(part_healthy)),clim1)
title('hypoactive - healthy')
colormap(s(3),CustomPastels)
colorbar
xlabel('TRs')
ylabel('ROIs')
%color plots - https://projects.susielu.com/viz-palette?colors=[%22#006593%22,%22#67bbba%22,%22#004a74%22,%22#000000%22,%22#fbc4ab%22,%22#f8ad9d%22,%22#f08080%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22normal%22
%plot the PC changing overtime (like a scatter plot)
figure
set(gcf,'color','w')
sig_avg_pc_hypo = squeeze(mean(part_hypo)).*sig_pc_hypo;
plot(mean(sig_avg_pc_hypo),'Color',[240 128 128]./255) %red color
sig_avg_pc_healthy = squeeze(mean(part_healthy)).*sig_pc_hypo;
hold on
plot(mean(sig_avg_pc_healthy),'Color',[0 101 147]./255) %blue color
xlabel('TRs')
ylabel('avg. PC across all ROIs')
legend('Hypoactive Del.','Healthy')

%plot the avg. PC (overtime) for each ROI on brain
figure
set(gcf,'color','w')
sig_avg_pc_hypo = squeeze(mean(part_hypo)).*sig_pc_hypo;
plot(mean(sig_avg_pc_hypo,3'),'Color',[240 128 128]./255) %red color
sig_avg_pc_healthy = squeeze(mean(part_healthy)).*sig_pc_hypo;
hold on
plot(mean(sig_avg_pc_healthy,3'),'Color',[0 101 147]./255) %blue color
xlabel('TRs')
ylabel('avg. PC across all ROIs')
legend('Hypoactive Del.','Healthy')

% variation in PC overtime

sig_std_pc_hypo = squeeze(std(part_hypo),[],3).*sig_pc_hypo);


%correlation of DRS score with participation coefficient values
reshape_all_part = reshape(part_all,69,502*210);

% correlation not working - as there are 'NaN' and 1's - creating
% dispreportionate numbers

[rho_drs_total,pval_drs_total]=corr(reshape_all_part,drs_total_calculated2,'rows','complete');
figure
set(gcf,'Color','w');
imagesc(reshape(rho_drs_total,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS')

figure
set(gcf,'Color','w');
subplot(1,4,1)
imagesc(reshape(rho_drs_total,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS All Subs')
[rho_drs_hypo,pval_drs_hypo]=corr(reshape_all_part(hypoactive_sub,:),drs_total_calculated2(hypoactive_sub,:),'rows','complete');
subplot(1,4,2)
imagesc(reshape(rho_drs_hypo,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS Hypo Subs')
[rho_drs_hyper,pval_drs_hyper]=corr(reshape_all_part(hyperactive_sub,:),drs_total_calculated2(hyperactive_sub,:),'rows','complete');
subplot(1,4,3)
imagesc(reshape(rho_drs_hyper,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS Hyper Subs')
[rho_drs_health,pval_drs_health]=corr(reshape_all_part(health_sub,:),drs_total_calculated2(health_sub,:),'rows','complete');
subplot(1,4,4)
imagesc(reshape(rho_drs_health,502,210));
xlabel('Time')
ylabel('Participation Coefficient')
title('Correlation between PC to total DRS Healthy Subs')


%% Cartographic profile
load('cartographic_all.mat')
cartographic_all([sub_remove],:,:,:)=[];
size(cartographic_all)
cart_hyper = cartographic_all(hyperactive_sub,:,:,:);
cart_hyper_mean = squeeze(mean(cart_hyper,4)); %avg. across time
cart_hypo = cartographic_all(hypoactive_sub,:,:,:);
cart_hypo_mean = squeeze(mean(cart_hypo,4));
cart_health = cartographic_all(health_sub,:,:,:);
cart_health_mean = squeeze(mean(cart_health,4));
cart_delirium = cartographic_all(delirium_sub,:,:,:);
cart_delirium_mean = squeeze(mean(cart_delirium,4));

figure
set(gcf,'color','w')
subplot(1,3,1)
imagesc(squeeze(mean(cart_health_mean)))
title('Avg. Cartograph Healthy')
xlabel('PC')
ylabel('Module degree')
subplot(1,3,2)
imagesc(squeeze(mean(cart_delirium_mean)))
title('Avg. Cartograph Delirium')
xlabel('PC')
ylabel('Module degree')
subplot(1,3,3)
diff_del_health = mean(cart_delirium_mean)-mean(cart_health_mean);
imagesc(squeeze(diff_del_health))
title('Avg. Cartograph Delirium - Health')
xlabel('PC')
ylabel('Module degree')

%sub group of del.
figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(squeeze(mean(cart_hypo_mean)))
title('Avg. Cartograph Hypoactive Del')
xlabel('PC')
ylabel('Module degree')
subplot(1,2,2)
imagesc(squeeze(mean(cart_hyper_mean)))
title('Avg. Cartograph Hyperactive Del')
xlabel('PC')
ylabel('Module degree')







%clustering the participation coeff. into network clusters
clustering_coef_bu();




% significant testing


%% Modularity
load('modularity_all.mat')
modularity_all([sub_remove],:)=[];
%1. Grps
mod_hyper = modularity_all(hyperactive_sub,:);
mod_hyper_mean = squeeze(mean(mod_hyper,1)); %avg. across time
mod_hypo = modularity_all(hypoactive_sub,:);
mod_hypo_mean = squeeze(mean(mod_hypo,1));
mod_health = modularity_all(health_sub,:);
mod_health_mean = squeeze(mean(mod_health,1));
mod_delirium = modularity_all(delirium_sub,:);
mod_delirium_mean = squeeze(mean(mod_delirium,1));


%stats
for i=1:size(modularity_all,2)
   [sig_mod_grp(i,:),pval_mod_grp(i,:)]=perm_code(mod_delirium(:,i),mod_health(:,i),1000);
end
%figures
figure
set(gcf,'color','w')
subplot(1,3,1)
imagesc(mod_delirium_mean')
title('Avg. Modularity Del.')
subplot(1,3,2)
imagesc(mod_health_mean')
title('Avg. Modularity Health')
subplot(1,3,3)
imagesc(sig_mod_grp)
title('Sig. Difference Time point')

%% 6. Modula Degree Z-score
load('modz_all.mat')
modz_all([sub_remove],:,:)=[];
%1. Grps
modz_hyper = modz_all(hyperactive_sub,:,:);
modz_hyper_mean = squeeze(mean(modz_hyper,3)); %avg. across time
modz_hypo = modz_all(hypoactive_sub,:,:);
modz_hypo_mean = squeeze(mean(modz_hypo,3));
modz_health = modz_all(health_sub,:,:);
modz_health_mean = squeeze(mean(modz_health,3));
modz_delirium = modz_all(delirium_sub,:,:);
modz_delirium_mean = squeeze(mean(modz_delirium,3));


%stats
for i=1:size(modz_all,2)
    for k=1:size(modz_all,3)
   [sig_modz_grp(i,k),pval_modz_grp(i,k)]=perm_code(modz_delirium(:,i,k),mod_health(:,i,k),1000);
    end
end

nROIs = size(modz_all,2);
for i=1:nROIs
   for k=1:size(modz_all,3)
   [sig_modz_grp(i,k),pval_modz_grp(i,k)]=perm_code(modz_delirium(:,i,k),modz_health(:,i,k),1000);
   end
end


%avg. PC between grps
figure
set(gcf,'Color','w')
subplot(1,2,1)
plot(mean(modularity_all(delirium_subs,:),1))
xlabel('Time')
ylabel('Modularity')
title('Avg. Modularity delirium subs')
subplot(1,2,2)
plot(mean(modularity_all(health_sub,:),1))
xlabel('Time')
ylabel('Modularity')
title('Avg. Modularity healthy subs')

% pearson correlation
avg_modularity = mean(modularity_all,2);
std_modularity = std(modularity_all,[],2);
[rho_avg_mod_drs_total,pval_avg_mod_drs_total]=corr(avg_modularity,drs_total_calculated2,'rows','complete');
[rho_std_mod_drs_total,pval_std_mod_drs_total]=corr(std_modularity,drs_total_calculated2,'rows','complete');
[rho_mod_drs_total,pval_mod_drs_total]=corr(modularity_all,drs_total_calculated2,'rows','complete');


figure
set(gcf,'Color','w');
plot(rho_mod_drs_total)
xlabel('Time')
ylabel('Rho')
title('Correlation between modularity to total DRS')
colormap(CustomColormap)

