function [loop_performance,results,h] = eta_loocv(gs_struct,groups,gs_min)

k_max = 1;
X = gs_struct.X;
groups = groups(:);

if nargin < 3
    gs_min = 3;
end

g_gs_idx = gs_struct.g_gs_idx;
keep_m = find(sum(g_gs_idx > 0, 2) >= gs_min);
g_gs_idx = g_gs_idx(keep_m,:);

[G,N] = size(X);
N_1 = sum(groups);
N_2 = N - N_1;

p_1 = 0.5; p_2 = 0.5;

% Leave-one-out cross validation

loop_performance = zeros(1,N);
h = zeros(1,N); h_k = zeros(1,N);

for fold = 1:N
    tic
    train_idx = true(1,N); train_idx(fold) = 0;
    
    F_train = X(:,train_idx); F_test = X(:,~train_idx);
    Y = groups(train_idx);
    N_fold = numel(Y);
    N_1_fold = sum(Y);
    N_2_fold = N_fold - N_1_fold;
    
%---- TRAIN ----%

[R_1,R_2,avg_var_1,avg_var_2,T_train] = rank_matching(F_train,g_gs_idx,'train',Y);
Delta = R_1 - R_2;

[R_1_new,R_2_new] = rank_matching(F_test,g_gs_idx,'test',T_train);
Delta_new = R_1_new - R_2_new;

Y_hat = double(Delta > 0);
Y_hat(Delta == 0) = 0.5;
eta = (sum(Y_hat(:,Y),2)/N_1_fold)*p_1+(sum(1-Y_hat(:,~Y),2)/N_2_fold)*p_2;
[eta,m] = sort(eta,'descend');

gamma_1 = sum(Delta(m,Y),2)/N_1_fold;
gamma_2 = sum(Delta(m,~Y),2)/N_2_fold;
rho = abs(gamma_1 - gamma_2);

theta = zeros(k_max,1); u = 1;
while u <= k_max
    eta_u = eta(1);
    tie = eta == eta_u;
    [rho_u,m_u] = sort(rho(tie),'descend');
    theta(u:u+numel(m_u)-1) = m(m_u);
    m(tie) = []; eta(tie) = []; rho(tie) = [];
    u = u+numel(m_u);
end
theta = theta(1:k_max);

h_u = Y_hat(theta(1),:);
eta_opt = (sum(h_u(:,Y),2)/N_1_fold)*p_1+(sum(1-h_u(:,~Y),2)/N_2_fold)*p_2;
k_opt = 1;
for k = 3:2:k_max
    h_u = mode(Y_hat(theta(1:k),:));
    eta = (sum(h_u(:,Y),2)/N_1_fold)*p_1+(sum(1-h_u(:,~Y),2)/N_2_fold)*p_2;
    if eta > eta_opt
        k_opt = k;
        eta_opt = eta;
    end
end
k = k_opt;
loop_performance(fold) = eta_opt;

%---- TEST ----%

Y_hat_new = double(Delta_new(theta) > 0);
Y_hat_new(Delta_new(theta) == 0) = 0.5;
h(fold) = Y_hat_new(1,:);
if k_max > 1
    h_k(fold) = mode(Y_hat_new(1:k));
end

toc
display([num2str(100*fold/N) '% Complete (' num2str(fold) ' folds)'])
end

TP = sum(h(groups));
FN = sum(1-h(groups));
FP = sum(h(~groups));
TN = sum(1-h(~groups));

accuracy = (TP/N_1)*p_1 + (TN/N_2)*p_2;
results(1,:) = [TP FN FP TN accuracy];

if k_max > 1
k_TP = sum(h_k(groups));
k_FN = sum(1-h_k(groups));
k_FP = sum(h_k(~groups));
k_TN = sum(1-h_k(~groups));

k_accuracy = (k_TP/N_1)*p_1 + (k_TN/N_2)*p_2;
results(2,:) = [k_TP k_FN k_FP k_TN k_accuracy];
end

