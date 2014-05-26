function [mu_R_1,mu_R_2,mu_diff_gs,mu_diff_stats] = mu_diff(gs_struct,groups,num_permut,gs_min)

X = gs_struct.X;
groups = groups(:);

if nargin < 4
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
[R_1,R_2] = rank_matching(X,g_gs_idx,'train',groups);

% Calculate rank conservation indices for class 1 and 2
mu_R_1 = mean(R_1(:,groups),2);
mu_R_2 = mean(R_2(:,~groups),2);

mu_diff = abs(mu_R_1-mu_R_2);

[mu_diff,sort_idx] = sort(mu_diff);
mu_diff_gs = flipud(gs(sort_idx)); 
mu_diff_G_gs = G_gs(sort_idx); mu_diff_G_pairs = G_pairs(sort_idx);
mu_R_1 = flipud(mu_R_1(sort_idx,:)); 
mu_R_2 = flipud(mu_R_2(sort_idx,:));

[mu_diff_u,first_idx] = unique(mu_diff,'first');
[mu_diff_u,last_idx] = unique(mu_diff,'last');
M_mu_diff_u = numel(mu_diff_u);

mu_diff_matrix = repmat(mu_diff,1,M_mu_diff_u);
mu_diff_u_matrix = repmat(mu_diff_u',M,1);
mu_diff_counts = sum(mu_diff_matrix < mu_diff_u_matrix)';
mu_diff_cdf = mu_diff_counts/M;

mu_diff_rand_cdf = zeros(M_mu_diff_u,1);
% converge_idx = [];
for p = 1:num_permut
    tic
    % First permute the class labels
    [Rint] = randperm(numel(groups)); 
    groups_rand = groups(Rint);
    
    [R_1_rand,R_2_rand] = rank_matching(X,g_gs_idx,'train',groups_rand);

    mu_R_1_rand = mean(R_1_rand(:,groups_rand),2);
    mu_R_2_rand = mean(R_2_rand(:,~groups_rand),2);
    
    mu_diff_rand = abs(mu_R_1_rand - mu_R_2_rand);
    
    mu_diff_rand_matrix = repmat(mu_diff_rand,1,M_mu_diff_u);
    mu_diff_rand_counts = sum(mu_diff_rand_matrix < mu_diff_u_matrix)';
    
    cdf_last = mu_diff_rand_cdf;
    mu_diff_rand_cdf = mu_diff_rand_cdf+mu_diff_rand_counts/(num_permut*M);
    
%     if ~sum(abs(mu_diff_rand_cdf - cdf_last)./cdf_last >= 0.01)
%         converge_idx = [converge_idx; p];
%         display(['Converges on ' num2str(p) 'th permutation.'])
%         break
%     end
    
    display([num2str(100*p/num_permut) '% Complete (' num2str(p) ' permutations)'])
    toc
end

% Compute P-values for rank conservation index values
mu_diff_pvals = zeros(M,1);
mu_diff_fdr = zeros(M,1);
for m = 1:M_mu_diff_u
    mu_diff_pvals(first_idx(m):last_idx(m)) = 1-mu_diff_rand_cdf(m);
    mu_diff_fdr(first_idx(m):last_idx(m)) = ...
        (1-mu_diff_rand_cdf(m))/(1-mu_diff_cdf(m));
end

mu_diff_stats = flipud([keep_gs(sort_idx),mu_diff_G_gs,mu_diff_G_pairs,...
    mu_diff,mu_diff_pvals,mu_diff_fdr]);
