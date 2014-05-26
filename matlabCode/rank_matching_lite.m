function [T,O,T_diff] = rank_matching_lite(X,g_gs_idx,groups)


M = size(g_gs_idx,1); % M is the total number of groups, or pathways

N_1 = sum(groups); 
N_2 = sum(~groups);

T = struct('T_1',{},'T_2',{}); % group orderings (rank templates)
O = struct('O_1',{},'O_2',{});
T_diff = zeros(M,1);

for m = 1:M
    g_m = g_gs_idx(m,g_gs_idx(m,:)>0);
    G_m = numel(g_m);
    G_m_choose_2 = G_m*(G_m-1)/2;
    
    T_0 = nchoosek(1:G_m,2);
    
    X_m = X(g_m,:); % Pull out rows in X corresponding to current group
    
    T_1 = (sum(X_m(T_0(:,1),groups) < X_m(T_0(:,2),groups),2)/N_1)>0.5;
    T_2 = (sum(X_m(T_0(:,1),~groups) < X_m(T_0(:,2),~groups),2)/N_2)>0.5;

    T(m).T_1 = T_1; T(m).T_2 = T_2;
    
    % Determine the actual order for each class
    for y = 1:2
        set = 1:G_m;
        T_0_tmp = T_0;
        T_tmp = eval(['T_',num2str(y)]);
        O_m = zeros(G_m,1);
        count = 1;
        while numel(set)
            i = set(count);
            i_lt = find(T_0_tmp(:,1) == i);
            T_i_lt = T_tmp(i_lt);
            i_gt = find(T_0_tmp(:,2) == i);
            T_i_gt = T_tmp(i_gt);
            gt_i = sum(T_i_lt) + sum(~T_i_gt);
            if gt_i == 0
                set = setdiff(set,i);
                O_m(G_m-numel(set)) = i;
                T_0_tmp([i_lt;i_gt],:) = [];
                T_tmp([i_lt;i_gt]) = [];
                count = 1;
            else count = count+1;
            end
        end
        eval(['O(m).O_',num2str(y),' = O_m;']);
    end
    
    T_diff(m) = sum(T_1~=T_2)/G_m_choose_2;
end



    