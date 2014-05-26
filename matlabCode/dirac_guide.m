
% ::: Load gene expression data :::

load prostate_GDS2545_m_nf

% Contains the following variables*:

% E_log10_QN - normalized gene expression matrix
% names - gene names
% groups - vector of class labels; 1 = class 1, 0 = class 2

% * see document "Variable Key" for more details on each

% Use the following commands to remove any rows containing missing values
if sum(sum(isnan(E_log10_QN)))
    notnan = find(sum(isnan(E_log10_QN),2)==0);
    names = names(notnan);
    E_log10_QN = E_log10_QN(notnan,:);
end

% If you wanted to create a specific expression matrix for each class, you
% could do the following:
E_MT = E_log10_QN(:,groups); % pulls out columns belonging to class 1
E_NP = E_log10_QN(:,~groups); % pulls out columns belonging to class 2

% In case you have separate expression matrices for two different classes
% but wanted to combine them into a typical binary dataset (for which most
% of the DIRAC functions are designed), you could use commands like these:
E_log10_QN = [E_NP,E_MT];
groups = [false(size(E_NP,2),1);true(size(E_MT,2),1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ::: Load pathway data :::

load gs_definitions biocarta_gs_defs

% Loads the matrix biocarta_gs_defs with genes grouped into pathways

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ::: Single phenotype rank conservation :::

% If you only wanted to find the rank conservation indices (mu_R) for one
% class at a time, you would follow these steps:

% Create the gene set data structure (gs_struct*) based on the expression
% data, the corresponding gene names, and the pathway data
gs_struct = gs_match_id(E_MT,names_new,biocarta_gs_defs);

% Calculate rank matching scores for all samples to the one class template
R_MT = rank_matching(gs_struct.X,gs_struct.g_gs_idx,'train',[]);

% Calculate rank conservation indices for each pathway
mu_R_MT = mean(R_MT,2);

% * see document "Variable Key" for more details on gs_struct components

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ::: Working with a typical binary dataset :::

% The following functions are designed to work for binary phenotype
% datasets (i.e., comparing two classes).  

% Typically when working with multiple datasets, I'll set up a script to
% load each, perform all calculations, and save the results to a structure.
% The value of d indicates the current dataset, and serves as an index for
% where to store data in the structure.

d = 1

% If you were working with multiple datasets, this is where you would load
% the data (as shown above). See any of the scripts ending in "_results"
% for an example.

% Create the gene set data structure (gs_struct)
gs_struct = gs_match_id(E_log10_QN,names_new,biocarta_gs_defs);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% NOTE: The "rank_matching" function is typically not used directly with a
% dataset, but it is at the core of DIRAC and is thus called by all other
% functions.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ::: Rank conservation indices :::

% Calculate rank conservation indices for class 1 and class 2 for all
% pathways. Pathway number, size in terms of both genes and pairs, and the
% rank conservation index are listed in the "_stats" variables.
% Corresponding pathway names are listed in "_gs" variables.
[mu_R_1_gs,mu_R_1_stats,mu_R_2_gs,mu_R_2_stats] =  mu_R(gs_struct,groups);

% If applicable, store results in the appropriate structure.
mu_R_struct(d).name = 'metastatic_normal'
mu_R_struct(d).mu_R_1_gs = mu_R_1_gs;
mu_R_struct(d).mu_R_1_stats = mu_R_1_stats;
mu_R_struct(d).mu_R_2_gs = mu_R_2_gs;
mu_R_struct(d).mu_R_2_stats = mu_R_2_stats;
save mu_R_struct_prostate mu_R_struct

% Determine most differentially regulated pathways; that is, those pathways
% with the greatest difference in rank conservation index between class 1
% and class 2. The "mu_diff" function uses random permutations to determine
% a P-value for each difference value.
[mu_R_1,mu_R_2,mu_diff_gs,mu_diff_stats] = mu_diff(gs_struct,groups,1);

% If applicable, store results in the appropriate structure.
mu_diff_struct(d).name = 'metastatic_normal'
mu_diff_struct(d).gs = mu_diff_gs;
mu_diff_struct(d).mu_R = [mu_R_1,mu_R_2];
mu_diff_struct(d).mu_diff = mu_diff_stats;
save mu_diff_struct_prostate mu_diff_struct

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ::: Rank difference scores :::

% Determine most differentially expressed pathways based on rank difference
% score and calculate associated statistics and accuracies
[eta_gs,eta_stats] = eta_fdr(gs_struct,groups,1);

% Evaluate classification performance in cross-validation.
[classifiers,results] = eta_loocv(gs_struct,groups,1);

% If applicable, store results in the appropriate structure.
eta_struct(d).name = 'metastatic_normal'
eta_struct(d).gs = eta_gs;
eta_struct(d).stats = eta_stats;
eta_struct(d).accuracy = results(5);
save eta_struct eta_struct
clear
clc

