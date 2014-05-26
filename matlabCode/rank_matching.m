function [R_1,R_2,avg_var_1,avg_var_2,T_train,Z_train] = rank_matching(X,g_gs_idx,varargin)

if numel(varargin)
    methods = {'train','test'};
    pname = varargin{1}; pval = varargin{2};
    method = find(strncmpi(pname,methods,numel(pname)));
    if isempty(method)
        error('Unknown method %s.',pname);
    end
    switch(method)
        case 1 % train
            if numel(pval)
                groups = pval;
            else groups = true(size(X,2),1);
            end
        case 2 % test
            T_1_2 = pval;
            groups = 1;
    end
    method = methods{method};
else
    error('Please specify method.');
end


M = size(g_gs_idx,1); % M is the total number of groups, or pathways
N = size(X,2); % N is samples

switch method
    case 'train'
    N_1 = sum(groups); 
    N_2 = sum(~groups);
end

R_1 = zeros(M,N); R_2 = zeros(M,N);
avg_var_1 = zeros(M,1); avg_var_2 = zeros(M,1);
T_train = struct('T_1',{},'T_2',{}); % group orderings (rank templates)
Z_train = struct('Z_1',{},'Z_2',{}); % (subnetwork rankings)
for m = 1:M
    g_m = g_gs_idx(m,g_gs_idx(m,:)>0);
    G_m = numel(g_m);
    
    G_m_choose_2 = G_m*(G_m-1)/2;
    T_0 = nchoosek(1:G_m,2);
    
    X_m = X(g_m,:); % Pull out rows in X corresponding to current group
    
    X_m_var_1 = var(X_m(:,groups)); X_m_var_2 = var(X_m(:,~groups));
    avg_var_1(m) = mean(X_m_var_1,2); avg_var_2(m) = mean(X_m_var_2,2);
    
    switch method
        case 'train'
            T_1 = (sum(X_m(T_0(:,1),groups) < X_m(T_0(:,2),groups),2)/N_1)>0.5;
            T_2 = (sum(X_m(T_0(:,1),~groups) < X_m(T_0(:,2),~groups),2)/N_2)>0.5;

            T_train(m).T_1 = T_1; T_train(m).T_2 = T_2;
         
        case 'test'
            T_1 = T_1_2(m).T_1;
            T_2 = T_1_2(m).T_2;
    end
    
    Z = X_m(T_0(:,1),:) < X_m(T_0(:,2),:);
    
    switch method
        case 'train'
            Z_train(m).Z_1 = Z(:,groups); Z_train(m).Z_2 = Z(:,~groups);
    end
    
    R_1(m,:) = sum(Z == repmat(T_1,1,N),1)/G_m_choose_2;
    R_2(m,:) = sum(Z == repmat(T_2,1,N),1)/G_m_choose_2;
end

switch method
    case 'train'
        if ~N_2
            R_2 = [];
        end
end