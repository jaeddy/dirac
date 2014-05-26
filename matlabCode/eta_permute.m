function [eta_gs,eta_stats] = eta_permute(gs_struct,groups,num_permut,gs_min)

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
[R_1,R_2,avg_var_1,avg_var_2,T_train] = ...
    rank_matching(X,g_gs_idx,'train',groups);

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

eta_rand_counts = zeros(M,1);

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
    Delta_rand = R_1_rand - R_2_rand;
    
    Y_hat_rand = double(Delta_rand > 0);
    Y_hat_rand(Delta_rand == 0) = 0.5;
    eta_rand = (sum(Y_hat_rand(:,groups_rand),2)/N_1)*p_1...
        +(sum(1-Y_hat_rand(:,~groups_rand),2)/N_2)*p_2;
    
    eta_rand_counts = eta_rand_counts ...
        + (eta_rand >= eta);
       
    display([num2str(100*p/num_permut) '% Complete (' num2str(p) ' permutations)'])
    toc
end
    
eta_pvals = eta_rand_counts./num_permut;

% Use the Storey method to estimate FDR and q-values
[eta_fdr,eta_qvals] = mafdr(eta_pvals);


[eta_pvals,sort_idx] = sort(eta_pvals,'descend');
eta = eta(sort_idx);
eta_gs = flipud(gs(sort_idx)); 
eta_G_gs = G_gs(sort_idx); eta_G_pairs = G_pairs(sort_idx);
eta_T_diff = T_diff(sort_idx);
% eta_pvals = eta_pvals(sort_idx);
eta_fdr = eta_fdr(sort_idx);
eta_qvals = eta_qvals(sort_idx);

eta_stats = flipud([keep_gs(sort_idx),eta_G_gs,eta_G_pairs,...
    eta_T_diff,eta,eta_pvals,eta_fdr,eta_qvals]);
