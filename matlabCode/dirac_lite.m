
function [results_table,T,O] = dirac_lite(X,names,groups,geneset_defs_file,gs_min)

% INPUTS:
% - X: gene expression matrix (rows = genes, columns = samples)
% - names: gene names
% - groups: vector of class labels; 1 = class 1, 0 = class 2
% - geneset_defs_file: filename of gene set definitions (e.g.,
%       gs_definitions.mat)
% - gs_min: minimum gene set size (e.g., 3)

% OUTPUTS:
% - results_tables: formatted results for all gene sets
% - T: pairwise rank templates for each class for each gene set
% - O: actual ordering of genes for each class for each gene set

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gene Set Mapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gs_struct = gs_match_id(X,names,geneset_defs_file);

%%%%%%%%%%%%%%%%%%%%%%%%%% Template difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the gene set data structure (gs_struct*) based on the expression
% data, the corresponding gene names, and the gene set data
g_gs_idx = gs_struct.g_gs_idx;
G_m = sum(g_gs_idx > 0, 2);
keep_gs = find(G_m >= gs_min);
G_m = G_m(keep_gs); G_m_choose_2 = G_m.*(G_m-1)./2;
g_gs_idx = g_gs_idx(keep_gs,:);
gs = gs_struct.gs(keep_gs);

% Determine rank templates for each sample and calculate Hamming distance
% between templates for each gene set
[T,O,T_diff] = rank_matching_lite(X,g_gs_idx,groups);

% Sort and organize results
[T_diff_s,s_idx] = sort(T_diff,'descend');
gs = gs(s_idx);
G_m = G_m(s_idx);
G_m_choose_2 = G_m_choose_2(s_idx);

T = T(s_idx)';
O = O(s_idx)';

results_table = [gs,strtrim(cellstr(num2str(G_m))),...
    strtrim(cellstr(num2str(G_m_choose_2))),...
    strtrim(cellstr(num2str(T_diff_s)))];
results_table = [{'Name','Num. genes','Num. pairs','Template diff'};...
    results_table];

