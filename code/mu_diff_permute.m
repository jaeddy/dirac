function [mu_R_1,mu_R_2,mu_diff_gs,mu_diff_stats] = mu_diff_permute(gs_struct,groups,num_permut,gs_min)

X = gs_struct.X;
groups = groups(:);
N = size(X,2);
N_1 = sum(groups);
N_2 = N-N_1;

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

mu_diff_rand_counts = zeros(M,1);

% For small sample sizes, search through all possible permutations (if
% fewer than user-defined number of permutations
if nchoosek(N,N_1) < num_permut
	fullPermute = 1;
	n = 1:N;
	permuteMat = nchoosek(n,N_1);
	num_permut = size(permuteMat,1);
else fullPermute = 0;
end

for p = 1:num_permut
    tic
    % First permute the class labels
    if fullPermute
    	Rint = permuteMat(p,:);
    	Rint = [Rint, setdiff(n,Rint)];
    else [Rint] = randperm(numel(groups));
    end
    groups_rand = groups(Rint);
    
    [R_1_rand,R_2_rand] = rank_matching(X,g_gs_idx,'train',groups_rand);

    mu_R_1_rand = mean(R_1_rand(:,groups_rand),2);
    mu_R_2_rand = mean(R_2_rand(:,~groups_rand),2);
    
    mu_diff_rand = abs(mu_R_1_rand - mu_R_2_rand);
    
    mu_diff_rand_counts = mu_diff_rand_counts ...
    	+ (mu_diff_rand >= mu_diff);
    
    display([num2str(100*p/num_permut) '% Complete (' num2str(p) ' permutations)'])
    toc
end

% Compute P-values for rank conservation index values
mu_diff_pvals = mu_diff_rand_counts./num_permut;

% Use the Storey method to estimate FDR and q-values
[mu_diff_fdr,mu_diff_qvals] = mafdr(mu_diff_pvals);

[mu_diff_pvals,sort_idx] = sort(mu_diff_pvals,'descend');
mu_diff = mu_diff(sort_idx);
mu_diff_gs = flipud(gs(sort_idx)); 
mu_diff_G_gs = G_gs(sort_idx); mu_diff_G_pairs = G_pairs(sort_idx);
% mu_diff_pvals = mu_diff_pvals(sort_idx);
mu_diff_fdr = mu_diff_fdr(sort_idx);
mu_diff_qvals = mu_diff_qvals(sort_idx);

mu_R_1 = flipud(mu_R_1(sort_idx));
mu_R_2 = flipud(mu_R_2(sort_idx));

mu_diff_stats = flipud([keep_gs(sort_idx),mu_diff_G_gs,mu_diff_G_pairs,...
    mu_diff,mu_diff_pvals,mu_diff_fdr,mu_diff_qvals]);
