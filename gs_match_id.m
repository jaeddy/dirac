function gs_struct = gs_match_id(X,g,geneset_defs_file)

geneset_defs_file(sum(~cellfun('isempty',geneset_defs_file(:,3:end)),2)>500,:) = [];
geneset_defs_file(sum(~cellfun('isempty',geneset_defs_file(:,3:end)),2)<1,:) = [];

gs = geneset_defs_file(:,1); 
gs_id = geneset_defs_file(:,2);
g_gs_file = geneset_defs_file(:,3:end); 
[M,G_m_max] = size(g_gs_file); 

[G,N] = size(X); 
[g_sorted,sort_idx] = sortrows(g);  % Sort gene names alphabetically.
X_sorted = X(sort_idx,:); % Update row order of gene expression matrix.

% Collapse gene expression matrix into rows of unique genes (i.e., exclude
% duplicate gene IDs); if duplicate IDs exist, use the max, per sample
% expression in the collapsed matrix.

[g_unique,first_idx] = unique(g_sorted,'first'); 
[g_unique,last_idx] = unique(g_sorted,'last'); 
G_unique = numel(g_unique); 

singleton = first_idx == last_idx; % Keep track of genes with only one instance.
multiple = [first_idx(~singleton) last_idx(~singleton)]; % Define the idx ranges for duplicates.
G_multiple = size(multiple,1); 

X_multiple = zeros(G_multiple,N); 
for i = 1:G_multiple 
    X_multiple(i,:) = max(X_sorted(multiple(i,1):multiple(i,2),:),[],1); 
end

X = zeros(G_unique,N); 
X(singleton,:) = X_sorted(first_idx(singleton),:); 
X(~singleton,:) = X_multiple; 

tic
g_gs_0 = repmat({''},M,2*G_m_max); 

for i = 1:M
    g_m = g_gs_file(i,~cellfun('isempty',g_gs_file(i,:)));
    multiple = ~cellfun('isempty',regexp(g_m,'\s')); 
    g_append = strcat(g_m(multiple),'///');
    if numel(g_append)
        g_multiple = strtrim(regexp(cell2mat(g_append),'///','split'));
    else g_multiple = {};
    end
    g_m = union(g_m(~multiple),g_multiple(~cellfun('isempty',g_multiple)));
    g_gs_0(i,1:numel(g_m)) = g_m; 
end
g_gs_0(:,sum(~cellfun('isempty',g_gs_0))==0) = []; 
g_gs_total = unique(g_gs_0); 

% Only search through gene names in the definitions file.
[g_intersect,int_idx] = intersect(g_unique,g_gs_total); 

g_gs = repmat({''},size(g_gs_0)); 
g_gs_idx = zeros(size(g_gs_0)); 
for i = 1:M
    g_m = g_gs_0(i,~cellfun('isempty',g_gs_0(i,:)));
    [g_m_intersect,int_idx_m] = intersect(g_intersect,g_m); 
    g_gs(i,1:numel(g_m_intersect)) = g_m_intersect; 
    g_gs_idx(i,1:numel(int_idx_m)) = int_idx(int_idx_m); 
end
toc

[g_gs_idx,unique_idx] = unique(g_gs_idx,'rows');
g_gs = g_gs(unique_idx,:); 
G_gs = sum((g_gs_idx > 0),2); 
G_gs_matched = numel(unique(g_gs)); 

g_gs_0 = g_gs_0(unique_idx,:);
G_gs_0 = sum(~cellfun('isempty',g_gs_0),2); 
g_gs_total = unique(g_gs_0);
G_gs_total = numel(g_gs_total); 
gs = gs(unique_idx); 
gs_id = gs_id(unique_idx);

g_gs_match_rate = [G_gs, G_gs_0, G_gs./G_gs_0]; 

gs_struct = struct('gs',{gs},...
    'gs_id',{gs_id},...
    'G_gs_total',{G_gs_total},...
    'G_gs_matched',{G_gs_matched},...
    'g_gs',{g_gs},...
    'g_gs_idx',{g_gs_idx},...
    'g_gs_match_rate',{g_gs_match_rate},...
    'X',{X},'g',{g_unique});

