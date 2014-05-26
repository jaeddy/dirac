function [mu_R_1_gs,mu_R_1_stats,mu_R_2_gs,mu_R_2_stats] = mu_R_fdr(gs_struct,groups,num_permut)

X = gs_struct.X;
groups = groups(:);
N = size(X,2);
N_1 = sum(groups);
N_2 = N-N_1;

p_1 = 0.5; p_2 = 0.5;

g_gs_idx = gs_struct.g_gs_idx;
keep_gs = find(sum(g_gs_idx > 0, 2) > 2);
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

% Sort conservation index values and corresponding gene set names/sizes
[mu_R_1,R_1_sort_idx] = sort(mu_R_1);
mu_R_1_gs = flipud(gs(R_1_sort_idx)); 
mu_R_1_G_gs = G_gs(R_1_sort_idx); mu_R_1_G_pairs = G_pairs(R_1_sort_idx);

[mu_R_2,R_2_sort_idx] = sort(mu_R_2);
mu_R_2_gs = flipud(gs(R_2_sort_idx)); 
mu_R_2_G_gs = G_gs(R_2_sort_idx); mu_R_2_G_pairs = G_pairs(R_2_sort_idx);

% Find unique conservation index values in class 1
[mu_R_1_u,first_idx_1] = unique(mu_R_1,'first');
[mu_R_1_u,last_idx_1] = unique(mu_R_1,'last');
M_mu_R_1_u = numel(mu_R_1_u);

% Find unique conservation index values in class 2
[mu_R_2_u,first_idx_2] = unique(mu_R_2,'first');
[mu_R_2_u,last_idx_2] = unique(mu_R_2,'last');
M_mu_R_2_u = numel(mu_R_2_u);

% Calculate probabilities for conservation indices in class 1
mu_R_1_matrix = repmat(mu_R_1,1,M_mu_R_1_u);
mu_R_1_u_matrix = repmat(mu_R_1_u',M,1);
mu_R_1_counts = sum(mu_R_1_matrix < mu_R_1_u_matrix)';
mu_R_1_cdf = mu_R_1_counts/M;

% Calculate probabilities for conservation indicies in class 2
mu_R_2_matrix = repmat(mu_R_2,1,M_mu_R_2_u);
mu_R_2_u_matrix = repmat(mu_R_2_u',M,1);
mu_R_2_counts = sum(mu_R_2_matrix < mu_R_2_u_matrix)';
mu_R_2_cdf = mu_R_2_counts/M;

mu_R_1_rand_cdf = zeros(M_mu_R_1_u,1);
mu_R_2_rand_cdf = zeros(M_mu_R_2_u,1);
converge_idx = [];
for p = 1:num_permut
    tic
    % First permute the class labels
    [Rint] = randperm(numel(groups)); 
    groups_rand = groups(Rint);
    
    [R_1_rand,R_2_rand] = rank_matching(X,g_gs_idx,'train',groups_rand);

    mu_R_1_rand = mean(R_1_rand(:,groups_rand),2);
    mu_R_2_rand = mean(R_2_rand(:,~groups_rand),2);
    
    mu_R_1_rand_matrix = repmat(mu_R_1_rand,1,M_mu_R_1_u);
    mu_R_1_rand_counts = sum(mu_R_1_rand_matrix < mu_R_1_u_matrix)';
    
    cdf_last = mu_R_1_rand_cdf;
    mu_R_1_rand_cdf = mu_R_1_rand_cdf+mu_R_1_rand_counts/(num_permut*M);
    
    if ~sum(abs(mu_R_1_rand_cdf - cdf_last)./cdf_last >= 0.01)
        converge_idx = [converge_idx; p];
        display(['Converges on ' num2str(p) 'th permutation.'])
    end
    
    mu_R_2_rand_matrix = repmat(mu_R_1_rand,1,M_mu_R_2_u);
    mu_R_2_rand_counts = sum(mu_R_2_rand_matrix < mu_R_2_u_matrix)';
    
    mu_R_2_rand_cdf = mu_R_2_rand_cdf+mu_R_2_rand_counts/(num_permut*M);
    
    display([num2str(100*p/num_permut) '% Complete (' num2str(p) ' permutations)'])
    toc
end

% Compute P-values for rank conservation index values
mu_R_1_pvals = zeros(M,1);
for m = 1:M_mu_R_1_u
    mu_R_1_pvals(first_idx_1(m):last_idx_1(m)) = 1-mu_R_1_rand_cdf(m);
end

mu_R_2_pvals = zeros(M,1);
for m = 1:M_mu_R_2_u
    mu_R_2_pvals(first_idx_2(m):last_idx_2(m)) = 1-mu_R_2_rand_cdf(m);
end

mu_R_1_stats = flipud(...
    [keep_gs(R_1_sort_idx),mu_R_1_G_gs,mu_R_1_G_pairs,mu_R_1,mu_R_1_pvals]);
mu_R_2_stats = flipud(...
    [keep_gs(R_2_sort_idx),mu_R_2_G_gs,mu_R_2_G_pairs,mu_R_2,mu_R_2_pvals]);


