function [eta_gs,eta_stats] = eta_fdr(gs_struct,groups,num_permut,gs_min)

X = gs_struct.X;
groups = groups(:);
N = size(X,2);
N_1 = sum(groups);
N_2 = N-N_1;

p_1 = 0.5; p_2 = 0.5;

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
[R_1,R_2,avg_var_1,avg_var_2,T_train] = rank_matching(X,g_gs_idx,'train',groups);

% Compare rank templates
T_diff = zeros(M,1);
for m = 1:M
    T_1 = double(T_train(m).T_1);
    T_2 = double(T_train(m).T_2);
    
    T_diff(m) = sum(T_2~=T_1)/numel(T_1);
end
% Calculate rank difference metric
Delta = R_1 - R_2;

Y_hat = double(Delta > 0);
Y_hat(Delta == 0) = 0.5;
eta = (sum(Y_hat(:,groups),2)/N_1)*p_1+(sum(1-Y_hat(:,~groups),2)/N_2)*p_2;

[eta,sort_idx] = sort(eta);
eta_gs = flipud(gs(sort_idx)); 
eta_G_gs = G_gs(sort_idx); eta_G_pairs = G_pairs(sort_idx);
eta_T_diff = T_diff(sort_idx);

[eta_u,first_idx] = unique(eta,'first');
[eta_u,last_idx] = unique(eta,'last');
M_eta_u = numel(eta_u);

eta_matrix = repmat(eta,1,M_eta_u);
eta_u_matrix = repmat(eta_u',M,1);
eta_counts = sum(eta_matrix < eta_u_matrix)';
eta_cdf = eta_counts/M;

eta_rand_cdf = zeros(M_eta_u,1);
% converge_idx = [];
for p = 1:num_permut
    tic
    % First permute the class labels
    [Rint] = randperm(numel(groups)); 
    groups_rand = groups(Rint);
    
    [R_1_rand,R_2_rand] = rank_matching(X,g_gs_idx,'train',groups_rand);
    Delta_rand = R_1_rand - R_2_rand;
    
    Y_hat_rand = double(Delta_rand > 0);
    Y_hat_rand(Delta_rand == 0) = 0.5;
    eta_rand = (sum(Y_hat_rand(:,groups_rand),2)/N_1)*p_1...
        +(sum(1-Y_hat_rand(:,~groups_rand),2)/N_2)*p_2;
    
    eta_rand_matrix = repmat(eta_rand,1,M_eta_u);
    eta_rand_counts = sum(eta_rand_matrix < eta_u_matrix)';
    
    cdf_last = eta_rand_cdf;
    eta_rand_cdf = eta_rand_cdf + eta_rand_counts/(num_permut*M);
    
%     if ~sum(abs(eta_rand_cdf - cdf_last)./cdf_last >= 0.01)
%         converge_idx = [converge_idx; p];
%         display(['Converges on ' num2str(p) 'th permutation.'])
%         break
%     end
    
    display([num2str(100*p/num_permut) '% Complete (' num2str(p) ' permutations)'])
    toc
end
    
eta_pvals = zeros(M,1);
eta_fdr = zeros(M,1);
for m = 1:M_eta_u
    eta_pvals(first_idx(m):last_idx(m)) = 1-eta_rand_cdf(m);
    eta_fdr(first_idx(m):last_idx(m)) = ...
        (1-eta_rand_cdf(m))/(1-eta_cdf(m));
end

eta_stats = flipud([keep_gs(sort_idx),eta_G_gs,eta_G_pairs,...
    eta_T_diff,eta,eta_pvals,eta_fdr]);
