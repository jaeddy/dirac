function [mu_R_1_gs,mu_R_1_stats,mu_R_2_gs,mu_R_2_stats] = mu_R(gs_struct,groups,gs_min)

X = gs_struct.X;
groups = groups(:);

if nargin < 3
    gs_min = 3;
end

g_gs_idx = gs_struct.g_gs_idx;
keep_gs = find(sum(g_gs_idx > 0, 2) >= gs_min);
g_gs_idx = g_gs_idx(keep_gs,:);
M = size(g_gs_idx,1);

gs = gs_struct.gs(keep_gs);
G_gs = gs_struct.g_gs_match_rate(keep_gs,1);
G_pairs = zeros(M,1);
for m = 1:M
    G_pairs(m) = G_gs(m)*(G_gs(m)-1)/2;
end

% Calculate rank matching scores for class 1 and 2
[R_1,R_2,avg_var_1,avg_var_2] = rank_matching(X,g_gs_idx,'train',groups);

% Calculate rank conservation indices for class 1 and 2
mu_R_1 = mean(R_1(:,groups),2);
mu_R_2 = mean(R_2(:,~groups),2);

% Sort conservation index values and corresponding gene set names/sizes
[mu_R_1,R_1_sort_idx] = sort(mu_R_1);
mu_R_1_gs = flipud(gs(R_1_sort_idx)); 
mu_R_1_G_gs = G_gs(R_1_sort_idx); mu_R_1_G_pairs = G_pairs(R_1_sort_idx);
avg_var_1 = avg_var_1(R_1_sort_idx);

[mu_R_2,R_2_sort_idx] = sort(mu_R_2);
mu_R_2_gs = flipud(gs(R_2_sort_idx)); 
mu_R_2_G_gs = G_gs(R_2_sort_idx); mu_R_2_G_pairs = G_pairs(R_2_sort_idx);
avg_var_2 = avg_var_2(R_2_sort_idx);

mu_R_1_stats = flipud(...
    [keep_gs(R_1_sort_idx),mu_R_1_G_gs,mu_R_1_G_pairs,avg_var_1,mu_R_1]);
mu_R_2_stats = flipud(...
    [keep_gs(R_2_sort_idx),mu_R_2_G_gs,mu_R_2_G_pairs,avg_var_2,mu_R_2]);


