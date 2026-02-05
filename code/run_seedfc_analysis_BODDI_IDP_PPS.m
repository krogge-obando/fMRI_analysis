%This analysis will conduct the analysis to see if there are any relationships
% between PPS, IDP and BBODI and the insula-whole brain analysis

%% Process behavioral spreadsheet
%% Load Excel behavioral table
behav_fname = '/data/neurogroup/parklab/social_data_fMRI_Oct25_BBODI.xlsx';
T = readtable(behav_fname, ...
              'Sheet', 'socialdata_Oct2025_clean', ...
              'VariableNamingRule', 'preserve');

%% Get all column labels (for reference)
table_labels = T.Properties.VariableNames;

%% Find rows where we want to keep subjects
rows_use = strcmp(T.("IfExclude_analysis_Oct25"), "No");
row_idx = find(rows_use);  % numeric version (optional)

%% Extract identifiers and covariates (keep same variable names)
sub_id      = T{:, 'ID'};       
MRI_id      = T{:, 'MRI_ID'};     
group       = T{:, 'Group'};       
sex         = T{:, 'Sex'};
hand        = T{:, 'Handness'};
age         = T{:, 'Age'};

%% Extract social_vr data (starting from column 12 to end)
social_vr_start_col = 11;
social_vr = T{:, social_vr_start_col:end};
table_labels_social_vr = table_labels(social_vr_start_col:end);



%% Now restrict to rows we’ll use in analysis
group     = group(rows_use);
sex       = sex(rows_use);
age       = age(rows_use);
hand      = hand(rows_use);
sub_id    = sub_id(rows_use);
MRI_id    = MRI_id(rows_use);
social_vr = social_vr(rows_use, :);


%% Convert group to boolean (1 = SZ, 0 = others)
group_bool = zeros(length(group), 1);
group_bool(strcmp(group, 'SZ')) = 1;

%% Display checks
fprintf('Included subjects: %d\n', sum(rows_use));
disp(table_labels);
disp(table_labels_social_vr);

%% Load Insula zscored maps into array

%load mni mask

maskfile = '/data/neurogroup/parklab/MNI152_T1_2mm_brain_Mask.nii';
mm = niftiread(maskfile);
bm = find(mm==1);


insula_maps_dir = '/data/neurogroup/parklab/PPS_study/ROI_maps/insula_2326'

maps_CO = containers.Map();
maps_SZ = containers.Map();

all_insula_mat = [];

% get maps for controls and SZ scan
for jj = 1: length(MRI_id)
    
    jj_MRI_id = MRI_id{jj};
    
    group_id = group{jj};
    
    nii_path    = [insula_maps_dir,'/',jj_MRI_id,'-scan01_INSULA_z.nii'];

    map_jj = niftiread(nii_path);
    
    all_insula_mat(jj,:) = map_jj(bm);
    
    clear map_jj

end



maps_SZ('insula') = all_insula_mat(group_bool==1,:);
maps_CO('insula') = all_insula_mat(group_bool==0,:);

% ---- Finish organizing maps ---- %
fprintf('Maps loaded and organized.\n');


%% Run t-tests between groups, for each seed
% Analysis to be done (for each ROI):
% 1. SZ and CO, each correlates with each of the social_vr measures
% 2. set group as a covariate, correlate all subjects' maps with each of the social_vr measures (SZ pre + CO only)
% 
seed_names = {'insula'};

% 
% configs for write NIFTI files
%derive masks, and rois
maskfile = '/data/neurogroup/parklab/MNI152_T1_2mm_brain_Mask.nii';
info = niftiinfo(maskfile);
info.Datatype = 'double';

stats_maps_dir = '/data/neurogroup/parklab/analysis/regressionFeb5_InsulaBODDI_IDP_PPS/';

mkdir(stats_maps_dir);


%% 1.1 CO separately
    % separately analyze controls only

beta_insula_con = {};

tmap_insula_con = {};
alpha_p = 0.05;
ifsaving_volume = true;


for cc=1:size(social_vr,2)
    cc_metric_name = table_labels_social_vr{cc};
    fprintf('Processing control group, social_vr measure: %s (%d of %d)\n', cc_metric_name, cc, size(social_vr,2));

    inds_con = find(group_bool==0);
    age_prep = age(inds_con);
    sex_prep = sex(inds_con);
    social_vr_prep = social_vr(inds_con,cc);

    % --- further filter the data, to exclude the nan values in social_vr
    not_nan_ind = ~isnan(social_vr_prep);
    num_valid = sum(not_nan_ind);
    social_vr_cc_prep = social_vr_prep(not_nan_ind);
% modeling covariates ind intercept
    XX = appendOnes([age_prep(not_nan_ind), ...
                 sex_prep(not_nan_ind), zscore(social_vr_cc_prep)]);
    
    n = size(XX,1);
    p = size(XX,2);
    c = [0 0 0 1]';
    tcritical = tinv(1-(alpha_p/2),n-p);

    tcritical_HC=tcritical;

    % Insula
    all_insula_con = maps_CO('insula');
    betas_tmp = pinv(XX)*all_insula_con(not_nan_ind,:);
    betas_insula = betas_tmp(end,:);
    err = all_insula_con(not_nan_ind,:) - XX*betas_tmp;
    errorVar = (1/(n-p))*sum(err.^2,1); % (n-p) degrees of freedom
    t_num = c'*betas_tmp; % numerator of t-statistic
    t_den = sqrt(errorVar*(c'*inv(XX'*XX)*c)); 
    t_stat_insula = t_num./t_den; % (1 x voxels) vector of t-statistics
    

    % make volumes
    vol_betas_insula = zeros(91,109,91);
    vol_betas_insula(bm) = betas_insula;


    % turn these into t-values and p-values (uncorrected) -> later run FDR
    vol_tmap_insula = zeros(91,109,91);
    vol_tmap_insula(bm) = t_stat_insula;



    if ifsaving_volume
        % Insula
        fname = [stats_maps_dir, cc_metric_name,'-beta_insula_CO.nii'];
        niftiwrite(vol_betas_insula,fname,info);
        fname = [stats_maps_dir, cc_metric_name,'-tmap_insula_CO.nii'];
        niftiwrite(vol_tmap_insula,fname,info);

    else
        beta_insula_con{cc} = vol_betas_insula;
        tmap_insula_con{cc} = vol_tmap_insula;

    end
    fprintf('--- Finished measure %s (Non-NaN subjects: %d)\n', cc_metric_name, num_valid);
    clear t_stat_insula;    
    clear betas_insula;
end

%% 1.2 SZ separately
    % separately analyze SZ only

beta_insula_sz = {};
tmap_insula_sz = {};
alpha_p = 0.05;
ifsaving_volume = true;


for cc=1:size(social_vr,2)
    cc_metric_name = table_labels_social_vr{cc};
    fprintf('Processing patient group, social_vr measure: %s (%d of %d)\n', cc_metric_name, cc, size(social_vr,2));

    inds_sz = find(group_bool==1);
    age_prep = age(inds_sz);
    sex_prep = sex(inds_sz);
    social_vr_prep = social_vr(inds_sz,cc);

    % --- further filter the data, to exclude the nan values in social_vr
    not_nan_ind = ~isnan(social_vr_prep);
    num_valid = sum(not_nan_ind);
    social_vr_cc_prep = social_vr_prep(not_nan_ind);

    XX = appendOnes([age_prep(not_nan_ind), ...
                 sex_prep(not_nan_ind), zscore(social_vr_cc_prep)]);
    
    n = size(XX,1);
    p = size(XX,2);
    c = [0 0 0 1]';
    tcritical = tinv(1-(alpha_p/2),n-p);

    tcritical_SP=tcritical;

    % insula
    all_insula_sz = maps_SZ('insula');
    betas_tmp = pinv(XX)*all_insula_sz(not_nan_ind,:);
    betas_insula = betas_tmp(end,:);
    err = all_insula_sz(not_nan_ind,:) - XX*betas_tmp;
    errorVar = (1/(n-p))*sum(err.^2,1); % (n-p) degrees of freedom
    t_num = c'*betas_tmp; % numerator of t-statistic
    t_den = sqrt(errorVar*(c'*inv(XX'*XX)*c)); 
    t_stat_insula = t_num./t_den; % (1 x voxels) vector of t-statistics
    
    % make volumes
    vol_betas_insula = zeros(91,109,91);
    vol_betas_insula(bm) = betas_insula;


    % turn these into t-values and p-values (uncorrected) -> later run FDR
    vol_tmap_insula = zeros(91,109,91);
    vol_tmap_insula(bm) = t_stat_insula;


    if ifsaving_volume
        % Insula
        fname = [stats_maps_dir, cc_metric_name,'-beta_insula_SZ.nii'];
        niftiwrite(vol_betas_insula,fname,info);
        fname = [stats_maps_dir, cc_metric_name,'-tmap_insula_SZ.nii'];
        niftiwrite(vol_tmap_insula,fname,info);
    else
        beta_insula_sz{cc} = vol_betas_insula;
        tmap_insula_sz{cc} = vol_tmap_insula;
    end
    fprintf('--- Finished measure %s (Non-NaN subjects: %d)\n', cc_metric_name, num_valid);
    clear t_stat_insula;    
    clear betas_insula;
end




%% 3. All subjects, with group as covariate
    % analyze all subjects, with group as a covariate

% combine pat/con (but covary for group status)
ifsaving_volume=true;
alpha_p = 0.05;

beta_insula_combined = {};
tmap_insula_combined = {};

% for cc=1:size(social_vr,2)
all_insula_maps_sz = maps_SZ('insula');
all_insula_maps_con = maps_CO('insula');
all_insula_mat = [all_insula_maps_sz; all_insula_maps_con];



for cc=1:size(social_vr,2) % inds in social_vr
    cc_metric_name = table_labels_social_vr{cc};
    fprintf('Processing control group, social_vr measure: %s (%d of %d)\n', cc_metric_name, cc, size(social_vr,2));


% --- further filter the data, to exclude the nan values in social_vr
    not_nan_ind = ~isnan(social_vr(:,cc));
    num_valid = sum(not_nan_ind);
    social_vr_cc_prep = social_vr(not_nan_ind,cc);
    XX = appendOnes([group_bool(not_nan_ind), age(not_nan_ind), sex(not_nan_ind), zscore(social_vr_cc_prep)]); % deleted hand for now because the social data is not completed
    
    n = size(XX,1);     
    p = size(XX,2);
    c = [0 0 0 0 1]';
%     tcritical = tinv(1-(0.001/2),n-p);
    tcritical = tinv(1-(alpha_p/2),n-p);
    tcritical_group=tcritical;

    % insula
    betas_tmp = pinv(XX)*all_insula_mat(not_nan_ind,:);
    betas_insula = betas_tmp(end,:);
    err = all_insula_mat(not_nan_ind,:) - XX*betas_tmp;
    errorVar = (1/(n-p))*sum(err.^2,1); % (n-p) degrees of freedom
    t_num = c'*betas_tmp; % numerator of t-statistic
    t_den = sqrt(errorVar*(c'*inv(XX'*XX)*c)); 
    t_stat_insula = t_num./t_den; % (1 x voxels) vector of t-statistics

    % make volumes
    vol_betas_insula = zeros(91,109,91);
    vol_betas_insula(bm) = betas_insula;

    % turn these into t-values and p-values (uncorrected) -> later run FDR
    vol_tmap_insula = zeros(91,109,91);
    vol_tmap_insula(bm) = t_stat_insula;

    if ifsaving_volume
        % insula
        fname = [stats_maps_dir, cc_metric_name,'-beta_insula_combined.nii'];
        niftiwrite(vol_betas_insula,fname,info);
        fname = [stats_maps_dir, cc_metric_name,'-tmap_insula_combined.nii'];
        niftiwrite(vol_tmap_insula,fname,info);
    else
        beta_insula_combined{cc} = vol_betas_insula;
        tmap_insula_combined{cc} = vol_tmap_insula;
    end
    fprintf('--- Finished measure %s (Non-NaN subjects: %d)\n', cc_metric_name, num_valid);
    clear t_stat_insula;    
    clear betas_insula;

end




