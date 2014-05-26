
gs_struct = gs_match_id(E_log10_QN,names,biocarta_gs_defs)

gs_struct = 
                 gs: {248x1 cell}
              gs_id: {248x1 cell}
         G_gs_total: 1296
       G_gs_matched: 1241
               g_gs: {248x87 cell}
           g_gs_idx: [248x87 double]
    g_gs_match_rate: [248x3 double]
                  X: [30744x68 double]
                  g: {30744x1 cell}
                  
X = gs_struct.X;

g_gs_idx = gs_struct.g_gs_idx;

g_m = 1:4
g_m =
     1     2     3     4
     
T_0 = nchoosek(1:numel(g_m),2)
T_0 =
     1     2
     1     3
     1     4
     2     3
     2     4
     3     4

mu_R_1 = mean(R_1(:,groups),2);
mu_R_2 = mean(R_2(:,~groups),2);
Delta = R_1 - R_2;

N_1 = sum(groups)
N_1 =
    37
    
N_2 = sum(~groups)
N_2 =
    31
    
eta = (sum(Delta(:,groups)>0,2)/N_1)*0.5 + (sum(Delta(:,~groups)<0,2)/N_2)*0.5;