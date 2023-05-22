% Script to define all Inputs for PLS for Medical Image Processing Toolbox
% 
% Parameters to be defined in this script are:
%
%   - input : struct containing input data for the analysis
%       - .X             : N x M matrix, N is #subjects, M is #imaging variables
%       - .Y             : N x B matrix, B is #behaviors
%       - .grouping  : N x 1 vector, subject grouping for PLS analysis
%                               e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                               subject 3 belongs to group 2.
%       - [.group_names]: Names of the groups (optional)
%       - [.behav_names]: Names of the behavior variables (optional) 
%   - pls_opts : options for the PLS analysis
%       - .nPerms              : number of permutations to run
%       - .nBootstraps         : number of bootstrapping samples to run
%       - .normalization_img   : normalization options for imaging data
%       - .normalization_behav : normalization options for behavior data
%              0 = no normalization
%              1 = zscore across all subjects
%              2 = zscore within groups (default)
%              3 = std normalization across subjects (no centering)
%              4 = std normalization within groups (no centering)
%       - .grouped_PLS : binary variable indicating if groups
%                                should be considered when computing R
%              0 = PLS will computed over all subjects
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices ( as in conventional behavior PLS)
%       - .grouped_perm : binary variable indicating if groups should be 
%                  considered during the permutations
%              0 = permutations ignoring grouping
%              1 = permutations within group
%       - .grouped_boot : binary variable indicating if groups should be 
%              considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - .boot_procrustes_mod : mode for bootstrapping procrustes transform
%              1 = standard (rotation computed only on U)
%              2 = average rotation of U and V
%       - .save_boot_resampling : binary variable indicating if bootstrap
%              resampling data should be saved or not
%              0 = no saving of bootstrapping resampling data
%              1 = save bootstrapping resampling data
%       - .behav_type        : Type of behavioral analysis
%              'behavior' for standard behavior PLS
%              'contrast' to simply compute contrast between two groups
%              'contrastBehav' to combine contrast and behavioral measures)
%              'contrastBehavInteract' to also consider group-by-behavior interaction effects
%   - save_opts: Options for result saving and plotting
%       - .output_path   : path where to save the results
%       - .prefix        : prefix of all results files (optional)
%       - .img_type      : Specify how to plot the results
%              'volume' for voxel-based data in nifti Format - results 
%                       will be displayed as bootstrap ratios in a brain map
%              'corrMat' for ROI-to-ROI correlation matrix - results will 
%                       be displayed as bootstrap ratios in a correlation matrix
%              'barPlot' for any type of brain data in already vectorized 
%                       form - results will be displayed as barplots
%       - .mask_file     : gray matter mask, only required if img_type='volume'
%       - .BSR_thres : 2x1 vector with negative and positive
%                          thresholds for bootstrap ratio map visualization,
%                          only required if imagingType='volume'
%       - .struct_file   : filename of structural file for background
%                          volume to overlay the results on, only required
%                          if imagingType='volume'
%       - .load_thres : 2x1 vector with negative and positive
%                          thresholds for loading map visualization,
%                          only required if imagingType='volume' 
%       - .grouped_plots : binary variable indicating if groups should be 
%                          considered during plotting
%              0 = plotting ignoring grouping
%              1 = plotting cosidering grouping
%       - .alpha         : significance level for LCs [default = 0.05]
%       - .plot_boot_samples : binary variable indicating if bootstrap
%                          samples should be plotted in bar plots
%       - .errorbar_mode : 'std' = plotting standard deviations
%                          'CI' = plotting 95% confidence intervals
%       - .hl_stable	 : binary variable indicating if stable bootstrap
%                          scores should be highlighted 


%% ---------- Input data ----------

%%%%%%%% LOAD YOUR DATA HERE %%%%%%%%
load('fc_fconn_Cos5_smth4_0exc_v02_wP070.mat');

pat_FC = cellstr(pat_ids);

load('Emo_behav_v02_0exc.mat');

X0 = fc_lined;

suf = '';

%X0 = cat(2, X0_occ, X0_dur);
%input.img_names = cat(2, names_occ, names_dur);

%[x, ia, ib] = intersect(pat_ids, pat_FC);
%%[x, ia, ib] = intersect(CovidCog, pat_FC);
%X0 = X0(ib, :);
%pat_FC = pat_FC(ib);

%[~, idx] = setdiff(pat_FC, CovidCog);
%X0(idx, :) = [];
%pat_FC(idx) = [];

%% EMO
Y0behav = [score_total, Fear, Disgust, Irritati1];
input.behav_names = {'Totale', 'Fear', 'Disgust', 'Irritation'};

%Y0behav = [score_total, Fear, Disgust, Irritati1, Age, Gender];
%input.behav_names = {'Totale', 'Fear', 'Disgust', 'Irritation', 'Age', ...
%                        'Gender'};

%Y0behav = [];
%input.behav_names = {};
%suf = '-memOnly';

%Y0behav = cat(2, Y0behav, [Age, Gender]);
%input.behav_names = cat(2, input.behav_names, {'Age', 'Gender'});
%suf = strcat('-wcov2', suf);
Y0behav = cat(2, Y0behav, [RLRI_DelayedTotalRecall]);
input.behav_names = cat(2, input.behav_names, {'Episodic memory'});
suf = strcat('-wmemo1', suf);

Y0behav = cat(2, Y0behav, [Age, NSC, Gender]);
input.behav_names = cat(2, input.behav_names, {'Age', 'NSC', 'Gender'});
suf = strcat('-wcov3', suf);

%% Handle missing values
is_nan = sum(isnan(Y0behav), 2) > 0;
%X0(is_nan, :) = [];
Y0behav(is_nan, :) = [];
pat_ids(is_nan) = [];

[x, ia, ib] = intersect(pat_ids, pat_FC);
X0 = X0(ib, :);
pat_FC = pat_FC(ib);

Y0behav = Y0behav(ia, :);
pat_ids = pat_ids(ia);
%CovidCog = CovidCog(ia);

%[~, idx] = setdiff(CovidCog, pat_FC);
%Y0behav(idx, :) = [];
%CovidCog(idx) = [];
%Gender(idx) = [];

%[x, ia, ib] = intersect(CovidCog, pat_FC);
%X0 = X0(ib, :);
%pat_FC = pat_FC(ib);

% Handle missing values
is_nan = sum(isnan(Y0behav), 2) > 0;
X0(is_nan, :) = [];
Y0behav(is_nan, :) = [];
pat_ids(is_nan) = [];

savename = [strcat('Emo-COVID', suf)];

GROUPING = 1;
grp_fltr = 0;
to_group = groups;

% Set to 1 for scatter plots (Warning: will require more than 1 min of run
% time and generate files larger than 2Gb)
SAVE_BOOTSTRAP = 0;

to_group = to_group(ia);
to_group(is_nan) = [];

%exclude = [];
% EMO: Remove FC from: [2, 10, 44] P077 -> 22
exclude = [22];

X0(exclude, :) = [];
pat_FC(exclude, :) = [];
Y0behav(exclude, :) = [];
to_group(exclude, :) = [];

group_names = {'Mild','Moderate', 'Severe'};
grp_names = {'_mild', '_moderate', '_severe'};

if grp_fltr
    X0(to_group ~= grp_fltr, :) = [];
    Y0behav(to_group ~= grp_fltr, :) = [];
    to_group(to_group ~= grp_fltr, :) = [];
    grp_suffix = char(grp_names(grp_fltr));
else
    grp_suffix = '';
end

% --- brain data ---
% Matrix X0 is typically a matrix with brain imaging data,
% of size subjects (rows) x imaging features (columns)
input.brain_data=X0;
%    '11','12', '13', '14', '15', '16', '17', '18', '19', '20', '21', ...
%    '22', '23', '24', '25', '26'};

% --- behavior data ---
% Matrix Y0 is a a matrix containing behavior data,
% of size subjects (rows) x behavior features (columns)
% Y0 will be constructed depending on the pls_opts.behav_type:
%   * if behav_type = 'contrast': behavData can be empty, Y0 will only depend
%     on the grouping information
%   * if behav_type = 'behavior': Y0 will be identical to behavData 
%   * if behav_type = 'contrastBehav', or 'contrastBehavInteract': 
%     Y0 will be constructed based on the grouping information and behavData
input.behav_data = double(Y0behav);

% --- grouping data ---
% subj_grouping: group assignment vector
% binary variable indicating the group, can contain multiple groups
input.grouping = to_group;

% --- Names of the groups ---
% here you can specify the names of the groups for the plots
if grp_fltr
    input.group_names = {grp_suffix};
else
    input.group_names=group_names;
end

%% ---------- Options for PLS ----------

% --- Permutations & Bootstrapping ---
pls_opts.nPerms = 1000; %initial is 200
pls_opts.nBootstraps = 500; %initial is 100
%pls_opts.nBootstraps = 380; %initial is 100

% --- Data normalization options ---
% 0: no normalization
% 1: zscore across all subjects
% 2: zscore within groups (default for grouped PLSC, see Krishnan et al.,2011)
% 3: std normalization across subjects (no centering)
% 4: std normalization within groups (no centering)
pls_opts.normalization_img = 1+GROUPING;
pls_opts.normalization_behav = 1+GROUPING;

% --- PLS grouping option ---
% 0: PLS will computed over all subjects
% 1: R will be constructed by concatenating group-wise covariance matrices
%     (as in conventional behavior PLS, see Krishnan et al., 2011)
pls_opts.grouped_PLS = GROUPING;

% --- Permutations grouping option ---
% 0: permutations ignoring grouping
% 1: permutations within group
pls_opts.grouped_perm = GROUPING;

% --- Bootstrapping grouping option ---
% 0: bootstrapping ignoring grouping
% 1: bootstrapping within group
pls_opts.grouped_boot = GROUPING;

% --- Mode for bootstrapping procrustes transform ---
% in some cases, rotation only depending on U results in extremely low
% standard errors and bootstrap ratios close to infinity
% in mode 2, we therefore compute the transformation matrix both on U and V
% 1: standard
% 2: average rotation of U and V
pls_opts.boot_procrustes_mod = 2;

% --- Save bootstrap resampling data? ---
% select whether bootstrapping resampling data should be saved (only
% recommended for few imaging dimensions)
pls_opts.save_boot_resampling=SAVE_BOOTSTRAP;

% --- Type of behavioral analysis ---
% 'behavior' for standard behavior PLS
% 'contrast' to simply compute contrast between two groups
% 'contrastBehav' to combine contrast and behavioral measures
% 'contrastBehavInteract' to also consider group-by-behavior interaction effects
pls_opts.behav_type = 'behavior';

%% ---------- Options for result saving and plotting ----------
% --- path where to save the results ---
path_to_results = '../../../outputs/PLS_results/';
save_opts.output_path = strcat(path_to_results, savename, grp_suffix, '_GRP'*GROUPING);

% --- prefix of all results files ---
% this is also the default prefix of the toolbox if you don't define
% anything
save_opts.prefix = sprintf('myPLS_%s_norm%d-%d',pls_opts.behav_type,...
    pls_opts.normalization_img, pls_opts.normalization_behav);

% --- Plotting grouping option ---
% 0: Plots ignoring grouping
% 1: Plots considering grouping
save_opts.grouped_plots = GROUPING;

% --- Significance level for latent components ---
save_opts.alpha = 0.1;
%save_opts.alpha = 0.05;

% --- Type of brain data ---
% Specify how to plot the results
% 'volume'  for voxel-based data in nifti Format - results will be 
%           displayed as bootstrap ratios in a brain map
% 'corrMat' for ROI-to-ROI correlation matrix - results will be displayed
%           as bootstrap ratios in a correlation matrix
% 'barPlot' for any type of brain data in already vectorized form - results 
%           will be displayed as barplots

% % uncomment the following to see example for volumetric plotting:
 %save_opts.img_type = 'volume'; 

% % uncomment the following to see example for correlation matrix plotting:
% input.brain_data=input.brain_data(:,1:30135); 
save_opts.img_type = 'corrMat';

% uncomment the following to see example for barplot figures:
%input.brain_data = input.brain_data(:,1:20); 
%save_opts.img_type = 'barPlot';
%save_opts.fig_pos_img = [440   606   560   192];

% --- Brain visualization thresholds ---
% (thresholds for bootstrap scores and loadings, only required if imagingType='volume' or imagingType='corrMat')
save_opts.BSR_thres = [-2.3 2.3]; % negative and positive threshold for visualization of bootstrap ratios
save_opts.load_thres = [-0.4 0.4]; % negative and positive threshold for visualization of loadings

% --- Brain mask ---
% (gray matter mask, only required if imagingType='volume')
save_opts.mask_file = 'example_mask.nii'; % filename of binary mask that will constrain analysis

% --- Structural template file for visualization ---
% (structural volume for background, only required if imagingType='volume')
save_opts.struct_file = 'example_struct.nii';

% --- Orientation of volumes in slice plots ---
% (only required if imagingType='volume')
save_opts.volume_orientation = 'axial'; %'axial','coronal','sagittal'

% --- Bar plot options ---
save_opts.plot_boot_samples = 1; % binary variable indicating if bootstrap samples should be plotted in bar plots
save_opts.errorbar_mode = 'CI'; % 'std' = plotting standard deviations; 'CI' = plotting 95% confidence intervals
save_opts.hl_stable = 1; % binary variable indicating if stable bootstrap scores should be highlighted

% --- Customized figure size for behavior bar plots ---
save_opts.fig_pos_behav = [440   606   320   192];
