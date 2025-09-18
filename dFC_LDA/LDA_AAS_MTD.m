%% Runninng LDA on the MTD (for the LC, nbM and PPN nuclei)

 % Load data
 load('AAS_MTD_permuted.mat')
% mtd_lc_avg: [nSubjects x nROIs x nTimeWindows]
data =mtd_ppn_avg;
name_suffix = 'mtd_lc';
[nSubj, nROIs, nTime] = size(data);
%% Option 2 - PCA on time-resolved data, X= subject * time, nROIs
%option 2 - % look at subject x time x roi (time componenet concate)
   % Each row = 1 time window from 1 subject
% Repeat label per time point
group_per_tp = repelem(bin_delirium_all, nTime); % [nSubj*nTime x 1], labels repeated per time

%export relevant data for python
X=reshape(data,nSubj*nTime,nROIs);
writematrix(X,'mtd_lc_avg_allsubstime.csv');
writematrix(X,'mtd_nbm_avg_allsubstime.csv');
writematrix(X,'mtd_ppn_avg_allsubstime.csv');
%writematrix(group_per_tp,'delirium_idx_repeattime.csv');

%save pc subject, nROIs x nTime
X2 = reshape(part_all,nSubj,nTime*nROIs);
writematrix(X2,'pc_flatten_roi.csv')

%% ------------- PCA
%pc_vec (coeff) - eigenvector, which is the 'loading' coefficient 
% pc_val (score)- pc scores, eigenvalues are simply the coefficients attached to eigenvectors, 
% which give the amount of variance carried in each Principal Component.
[pc_vec, pc_val, latent, ~, explained,mu] = pca(X); %mu, estimated mean of each variable in X

% Define normalized RGB colors
control_color = [253 205 154] / 255;   % Peach
delirium_color = [174 216 230] / 255;  % Light blue

%reshape first 2 PCs
pc_val_reshaped = reshape(pc_val(:, 1:2), [nSubj, nTime, 2]);
% reshape to subjects × time
pc1 = squeeze(pc_val_reshaped(:,:,1));  % [nSubj x nTime]
pc2 = squeeze(pc_val_reshaped(:,:,2));  % [nSubj x nTime]
pc1_mean = mean(pc1, 2);  % [nSubj x 1]
pc2_mean = mean(pc2, 2);  % [nSubj x 1]
% Group-wise average over time
nondel_pc1 = mean(pc_val_reshaped(bin_delirium_all==0,:,1), 1);
delirium_pc1 = mean(pc_val_reshaped(bin_delirium_all==1,:,1), 1);

% Plot
figure;
set(gcf, 'Color', 'w');
hold on;
scatter(pc1_mean(health_sub), pc2_mean(health_sub), 60, ...
    'MarkerFaceColor', control_color, 'MarkerEdgeColor', 'k');
scatter(pc1_mean(delirium_sub), pc2_mean(delirium_sub), 60, ...
    'MarkerFaceColor', delirium_color, 'MarkerEdgeColor', 'k');
xlabel(['PC 1 (' num2str(round(explained(1),1)) '%)']);
ylabel(['PC 2 (' num2str(round(explained(2),1)) '%)']);
legend({'Control','Delirium'}, 'Location', 'best');
title('Group Separation in PC Space');
axis square;
box on;
filename = sprintf('2_group_separation_pcs_%s.eps', name_suffix);
print('-depsc2', filename);  % '-depsc2' is color EPS

%% screee plot - percentage explained variance + cumulative variance
% Find where cumulative variance crosses 80%
cumVar = cumsum(explained);
idx80 = find(cumVar >= 80, 1, 'first'); %this is the number of PCs to 80%
var80 = cumVar(idx80);
% Plot
figure
set(gcf, 'Color', 'w')
yyaxis left
plot(explained, 'o-', 'LineWidth', 2);
ylabel('Individual Variance (%)');
yyaxis right
plot(cumVar, 's--', 'LineWidth', 2);
hold on;
plot(idx80, var80, 'rs', 'MarkerFaceColor', 'r');  % Red square at 80%
xline(idx80, '--r', 'LineWidth', 1);               % Vertical line at 80%
% Annotate
text(idx80 + 0.5, var80 - 5, sprintf('%.1f%% @ PC %d', var80, idx80), ...
     'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Cumulative Variance (%)');
xlabel('Principal Component');
title('PCA Variance Explained');
legend({'Individual', 'Cumulative', '80% threshold'}, 'Location', 'best');
grid on;
filename = sprintf('2_screeplot_cumulativevar_%s.eps', name_suffix);
print('-depsc2', filename);% '-depsc2' is color EPS



%% LDAs - on the time-features component
n_pcs = 30;
% define data classes (divisions - based upon the deliirum or non-delirious)
X1 = pc_val((group_per_tp==1),1:n_pcs); %delirium
X2 = pc_val(group_per_tp==0,1:n_pcs); %non-delirious
%Class sizes
N1 = size(X1,1);
N2 = size(X2,1);
%Class mean
Mu1 = mean(X1,1)';
Mu2 = mean(X2,1)';
%Data mean
Mu = (Mu1 + Mu2)./2;
% Between-class scatter martrix
Sb = N1.*(Mu1 - Mu)*(Mu1 - Mu)' + N2.*(Mu2 - Mu)*(Mu2 - Mu)';
% Within-class scatter matrices
S1 = cov(X1);
S2 = cov(X2);
Sw = S1 + S2; %aggregate within-class scatter
inv_Sw = inv(Sw);

G = inv_Sw*Sw; % Check inverse - QC for identity matrix

[eig_vec, eig_val] = eig(inv_Sw*Sb); % Eigendecomposition

D = real(diag(eig_val)); % QC

[~, eig_order] = sort(D,'descend');
lda_eig = eig_vec(:,eig_order); % Sort eigenvalues in descending order
% ---------------- Orthonormalize
% Not neccessary but but can be useful for eigenvector interpretation and
% data reconstruction
A = lda_eig;
% Modified Gram-Schmidt. [Q,R] = mgs(X);
% G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
[n,p] = size(A);
Q = zeros(n,p);
R = zeros(p,p);
Q(:,1) = A(:,1);
for kk = 2:p
    Q(:,kk) = A(:,kk);
    for ii = 1:kk-1
      R(ii,kk) = Q(:,ii)'*Q(:,kk);
      Q(:,kk) = Q(:,kk) - R(ii,kk)*Q(:,ii);
    end
    R(kk,kk) = norm(Q(:,kk))';
    Q(:,kk) = Q(:,kk)/R(kk,kk);
end
lda_eig = Q;
% LDA vectors in original space - loading of pc coefficient (pc
% eigenvector) * lda coefficient
lda_orig_space = pc_vec(:,1:n_pcs)*lda_eig;
% Project data into LDA subspace
x_lda_proj = pc_val(:,1:n_pcs)*lda_eig;



RB_surf_schaef_color(lda_orig_space(1:400,1),'~',[196 64 48],[12 125 121])


lda_1d = x_lda_proj(:,1);  % First LDA axis
group = bin_delirium_all;
lda_scores_reshaped = reshape(x_lda_proj, [nSubj,nTime,n_pcs]); 

% Split by group
lda_delirium = lda_scores_reshaped(:, delirium_sub, :);   % [210 x N_del x 8]
lda_nondel   = lda_scores_reshaped(:, health_sub, :);     % [210 x N_ctl x 8]

% Compute means over subjects (2nd dimension)
mean_delirium = mean(lda_delirium, 2);   % [210 x 1 x 8]
mean_control  = mean(lda_nondel, 2);     % [210 x 1 x 8]

% Compute standard error
stderr_delirium = std(lda_delirium, 0, 2) / sqrt(sum(delirium_sub)); % [210 x 1 x 8]
stderr_control  = std(lda_nondel, 0, 2) / sqrt(sum(health_sub));     % [210 x 1 x 8]

% Pull out first LDA axis (1st PC → LDA axis 1)
lda_dim = 1;  % or whichever LDA component you're interested in
m_del = squeeze(mean_delirium(:, 1, lda_dim));   % [210 x 1]
m_con = squeeze(mean_control(:, 1, lda_dim));    % [210 x 1]
se_del = squeeze(stderr_delirium(:, 1, lda_dim));  % [210 x 1]
se_con = squeeze(stderr_control(:, 1, lda_dim));   % [210 x 1]

% Plot
time = 1:nTime;
figure;
set(gcf, 'Color', 'w');
hold on;

% Plot shaded error bars
shadedErrorBar(m_del, se_del, ...
    'lineprops', {'-', 'Color', delirium_color, 'LineWidth', 2}); % Light blue

shadedErrorBar(m_con, se_con, ...
    'lineprops', {'-', 'Color', control_color, 'LineWidth', 2}); % Peach

xlabel('Time');
ylabel('LDA Score');
legend({'Delirium','Control'});
title('Mean LDA Score Over Time by Group - First LDA Axis');
grid on;


figure; 
set(gcf, 'Color', 'w'); 
hold on;

plot(time, mean_delirium(:,1,1), 'Color', [174 216 230]/255, 'LineWidth', 2);
plot(time, mean_control(:,1,1),  'Color', [253 205 154]/255, 'LineWidth', 2);

xlabel('Time');
ylabel('LDA Score');
legend({'Delirium','Control'});
title('Mean LDA Score Over Time by Group - based on 1st PC feature');
grid on;

%plotting the results onto a histogram
figure
histogram(x_lda_proj(delirium_sub,1)) %delirium sub
hold on
histogram(x_lda_proj(health_sub,1)) %health sub
title("histogram of the data projected into LDA space")






%% saving outputs
save('LDA_PCA_MTD.mat', 'lda_orig_space', 'lda_eig', 'pc_val','pc_vec', 'X','X1','X2','x_lda_proj')





%lda_orig_space - is 105420 x 10
nROIs = size(mtd_lc_avg,2);
nWindows = size(mtd_lc_avg,3);
lda_component_1 = reshape(lda_orig_space(:,1), [nROIs, nWindows]);  % [502 × 201]

roi_importance = mean(lda_component_1, 2);  % Can be positive or negative



lda_scores_matrix = reshape(x_lda_proj, nTime, nSubj,10);
%plot first LDA projection
RB_surf_schaef_color(lda_orig_space(1:400,1),'~',[196 64 48],[12 125 121])
figure;
hold on;
plot(mean(lda_scores_matrix(:,  health_sub,1), 2), 'b', 'LineWidth', 2); % Control
plot(mean(lda_scores_matrix(:, delirium_sub,1), 2), 'r', 'LineWidth', 2); % Delirium
xlabel('Time Window');
ylabel('LDA Score');
legend({'Control','Delirium'});
title('Temporal LDA Discriminant Trajectories');

figure;
bar(lda_orig_space(:,1));  % first LDA component
xlabel('ROI index');
ylabel('LDA weight');
title('ROI contributions to group separation (1st LDA component)');
grid on;


group_labels = repelem(bin_delirium_all, nTime);  % if scores are per time window

% Or average per subject if your scores are per subject
% lda_scores = mean(x_lda_proj(:,1), 2);

mean_delirium = mean(x_lda_proj(group_labels == 1,:));
mean_control  = mean(x_lda_proj(group_labels == 0,:));




