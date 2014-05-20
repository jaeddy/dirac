
function [results_tables,results_raw] = dirac(X,names,groups,...
    geneset_defs_file,geneset_defs_opt,gs_min,top_gs_M,num_permutations)

% INPUTS*:
% - X: gene expression matrix (rows = genes, columns = samples)
% - names: gene names
% - groups: vector of class labels; 1 = class 1, 0 = class 2
% - geneset_defs_file: filename of gene set definitions (e.g.,
%       gs_definitions.mat)
% - geneset_defs_opt: variable name of specific option for gene set
%       definitions (e.g., 'biocarta_gs_defs')
% - gs_min: minimum gene set size (e.g., 3)
% - top_gs_M: number of top gene sets to report
% - num_permutations: number of random permutations for calculating
%       p-values and false discovery rates (e.g., 1000)

% OUTPUTS:
% - results_tables: formatted results for top M gene sets
% - results_raw: all results stored in structures

% * see document "Variable Key" for more details on each

% EXAMPLE USAGE (or see "results_script.m"):
% - for provided metastatic vs. normal prostate dataset
%
% load('prostate_GDS2545_m_nf.mat')
% X = E;
% names = names_new;
% geneset_defs_file = 'gs_definitions.mat';
% geneset_defs_opt = 'biocarta_gs_defs';
% gs_min = 3;
% top_gs_M = 20;
% num_permutations = 10;

%[results_tables,results_raw] = dirac(X,names,groups,...
%    geneset_defs_file,geneset_defs_opt,gs_min,top_gs_M,num_permutations)


%%%%%%%%%%%%%%%%%%%%%%%%% Collect Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%

display(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gene Set Mapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(geneset_defs_file,geneset_defs_opt)
eval(['gs_defs = ', geneset_defs_opt, ';']);

display('Mapping gene IDs onto gene set defintions...')
gs_struct = gs_match_id(X,names,gs_defs);
results_raw.gs_struct = gs_struct;

display(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rank Matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the gene set data structure (gs_struct*) based on the expression
% data, the corresponding gene names, and the gene set data
g_gs_idx = gs_struct.g_gs_idx;
G_m = sum(g_gs_idx > 0, 2);
keep_gs = find(G_m >= gs_min);
G_m = G_m(keep_gs);
g_gs_idx = g_gs_idx(keep_gs,:);
results_raw.G_m = G_m;

if ~numel(top_gs_M)
    top_gs_M = numel(keep_gs);
end

% Calculate rank conservation indices for class 1 and class 2 for all
% gene sets.
[R_1,R_2,avg_var_1,avg_var_2,T_train] = ...
    rank_matching(X,g_gs_idx,'train',groups);
results_raw.T_train = T_train;

% * see document "Variable Key" for more details on gs_struct components

%%
%%%%%%%%%%%%%%%%%%%%%%% Tightly Regulated Networks %%%%%%%%%%%%%%%%%%%%%%%%

% Gene set number, size in terms of both genes and pairs, and the
% rank conservation index are listed in the "_stats" variables.
% Corresponding gene set names are listed in "_gs" variables.
[mu_R_1_gs,mu_R_1_stats,mu_R_2_gs,mu_R_2_stats] = ...
    mu_R(gs_struct,groups,gs_min);
results_raw.mu_R_1_gs = mu_R_1_gs;
results_raw.mu_R_1_stats = mu_R_1_stats;
results_raw.mu_R_2_gs = mu_R_2_gs;
results_raw.mu_R_2_stats = mu_R_2_stats;

col_headers = {'Network name','Num. genes','Num. pairs','Avg. variance',...
    'mu_R'};

num_col = size(mu_R_1_stats,2);
mu_R_1_stats_table = repmat({''},top_gs_M,num_col-1);
mu_R_2_stats_table = repmat({''},top_gs_M,num_col-1);
for i = 1:num_col-1
    mu_R_1_stats_table(:,i) = ...
        cellstr(num2str(mu_R_1_stats(1:top_gs_M,i+1)));
    mu_R_2_stats_table(:,i) = ...
        cellstr(num2str(mu_R_2_stats(1:top_gs_M,i+1)));
end

mu_R_1_table = [col_headers;[mu_R_1_gs(1:top_gs_M),mu_R_1_stats_table]];
mu_R_2_table = [col_headers;[mu_R_2_gs(1:top_gs_M),mu_R_2_stats_table]];

results_tables.Tightly_Regulated_Networks_Class_1 = mu_R_1_table;
results_tables.Tightly_Regulated_Networks_Class_2 = mu_R_2_table;

%%
%%%%%%%%%%%%%%%%%%%% Differentially Regulated Networks %%%%%%%%%%%%%%%%%%%%

% Determine most differentially regulated pathways; i.e., those gene sets
% with the greatest difference in rank conservation index between class 1
% and class 2. The "mu_diff" function uses random permutations to determine
% a P-value for each difference value.
[mu_R_1,mu_R_2,mu_diff_gs,mu_diff_stats] = ...
    mu_diff_permute(gs_struct,groups,num_permutations,gs_min);
results_raw.mu_R_1 = mu_R_1;
results_raw.mu_R_2 = mu_R_2;
results_raw.mu_diff_gs = mu_diff_gs;
results_raw.mu_diff_stats = mu_diff_stats;

num_sig = sum(mu_diff_stats(:,5)<.05);
num_sig = max(num_sig,top_gs_M);

col_headers = {'Network name','Num. genes','Num. pairs','mu_R in 1',...
    'mu_R in 2','Abs. difference in mu_R','P-value','FDR','Q-value'};

num_col = size(mu_diff_stats,2);
mu_diff_stats_table = repmat({''},num_sig,num_col-1);
for i = 1:num_col-1
    mu_diff_stats_table(:,i) = ...
        cellstr(num2str(mu_diff_stats(1:num_sig,i+1)));
end

mu_diff_stats_table = [mu_diff_stats_table(:,1:2),...
    cellstr(num2str(mu_R_1(1:num_sig))),...
    cellstr(num2str(mu_R_2(1:num_sig))),...
    mu_diff_stats_table(:,3:6)];

mu_diff_table = [col_headers;[mu_diff_gs(1:num_sig),mu_diff_stats_table]];

results_tables.Differentially_Regulated_Networks = mu_diff_table;
display(' ')

%%
%%%%%%%%%%%%%%%%%%%%%% Variably Expressed Networks %%%%%%%%%%%%%%%%%%%%%%%%

display('Generating null distribution of classification accuracies...')

% Determine most differentially expressed pathways based on rank difference
% score and calculate associated statistics and accuracies

% Generate null distribution of accuracies with random permutations
[eta_gs,eta_stats] = eta_permute(gs_struct,groups,num_permutations,gs_min);
results_raw.eta_gs = eta_gs;
results_raw.eta_stats = eta_stats;

num_sig = sum(eta_stats(:,6)<.05);
num_sig = max(num_sig,top_gs_M);

col_headers = {'Network name','Num. genes','Num. pairs',...
    'Template difference','Apparent accuracy','P-value','FDR','Q-value'};

num_col = size(eta_stats,2);
eta_stats_table = repmat({''},num_sig,num_col-1);
for i = 1:num_col-1
    eta_stats_table(:,i) = cellstr(num2str(eta_stats(1:num_sig,i+1)));
end

eta_table = [col_headers;[eta_gs(1:num_sig),eta_stats_table]];

results_tables.Variably_Expressed_Networks = eta_table;
display(' ')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Cross Validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display(['Performing cross validation to estimate classifation'...
    ' accuracy on future cases...']);

% Evaluate classification performance in cross-validation.
[training_accuracy,results,h] = eta_loocv(gs_struct,groups,gs_min);
results_raw.cv_results = results;
results_raw.predicted_classes = h';

col_headers = {'TP','FN','FP','TN','Accuracy'};
cv_table = [col_headers;strtrim(cellstr(num2str(results'))')];

results_tables.Classification_CV = cv_table;
display(' ')

display('Complete.')
