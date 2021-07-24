function [G_mn_R_p,G_mk_R_m,G_mo_R_in,G_mq_p,G_mr_m]=Compute_G_matrices(m_set,n1_set,n1_set_ind,n2_set,n2_set_ind,alpha_m,a_diff,beta1_m,beta0,beta2_m,L_mn_beta1,L_mk_beta2,L_mo_beta1,eig1_p,F_mq_p,eig2_m,F_mr_m,StrucParam)

a_diff_col = a_diff(StrucParam.N_Tr:1:2*StrucParam.N_Tr-1);
a_diff_row = a_diff(StrucParam.N_Tr:-1:1);
a_matr = toeplitz(a_diff_col, a_diff_row); % toeplitz matrix a_diff
a_sq = eye(StrucParam.N_Tr)+a_matr*a_matr; % 1+a_dif*a_dif

%% G_mn_R_p
if ~isempty(n1_set_ind)
    G_mn_R_p = zeros(length(m_set), length(n1_set));
    
    for m=1:length(m_set)
        for n=1:length(n1_set)
            
            G_mn_s = zeros(length(m_set),1); % preallocation
            
            for s=1:length(m_set)
                G_mn_s(s) = (a_matr(m,s) * alpha_m(s) - a_sq(m,s) * beta1_m(n1_set_ind(n))) * L_mn_beta1(s,n);
            end
            
            G_mn_R_p(m,n) = (StrucParam.Z0/(2*pi/StrucParam.lambda)/StrucParam.eps1) * sum(G_mn_s);
            
        end
    end
    
else G_mn_R_p = [];
end

%% G_mk_R_m
if ~isempty(n2_set_ind)
    G_mk_R_m = zeros(length(m_set),length(n2_set));
    
    for m=1:length(m_set)
        for k=1:length(n2_set)
            
            G_mk_s = zeros(length(m_set),1); % preallocation
            
            for s=1:length(m_set)
                G_mk_s(s) = (a_matr(m,s) * alpha_m(s) - a_sq(m,s) * beta2_m(n2_set_ind(k))) * L_mk_beta2(s,k);
            end
            
            G_mk_R_m(m,k) = (StrucParam.Z0/(2 * pi/StrucParam.lambda)/StrucParam.eps2) * sum(G_mk_s);
            
        end
    end
    
else G_mk_R_m = [];
end

%% G_mo_R_in
G_mo_R_in = zeros(length(m_set),1);
for m=1:length(m_set)
    
    G_mo_s = zeros(length(m_set),1); % preallocation
    
    for s=1:length(m_set)
        G_mo_s(s) = (a_matr(m,s) * alpha_m(s) - a_sq(m,s) * beta0) * L_mo_beta1(s);
    end
    
    G_mo_R_in(m,1) = (StrucParam.Z0/(2*pi/StrucParam.lambda)/StrucParam.eps1) * sum(G_mo_s);
end

%% G_mq_p
if length(n1_set) < StrucParam.N_Tr %condition that complex modes exist
    G_mq_p = zeros(length(m_set),StrucParam.N_Tr-length(n1_set));
    
    for m=1:length(m_set)
        for q=1:1:StrucParam.N_Tr-length(n1_set)
            
            G_mq_s = zeros(length(m_set),1); % preallocation
            
            for s=1:length(m_set);
                G_mq_s(s) = (a_matr(m,s) * alpha_m(s) - a_sq(m,s) * eig1_p(length(n1_set)+q)) * F_mq_p(s,q);
            end
            
            G_mq_p(m,q) = (StrucParam.Z0/(2*pi/StrucParam.lambda)/StrucParam.eps1) * sum(G_mq_s);
            
        end
    end
    
else G_mq_p = [];
end

%% G_mr_m
%{
a=length(n2_set)
b=length(m_set)
c=StrucParam.N_Tr-length(n2_set)
d=size(a_matr)
e=size(alpha_m)
f=size(a_sq)
g=size(eig2_m)
h=size(F_mr_m)
%}

if length(n2_set)<StrucParam.N_Tr
    G_mr_m = zeros(length(m_set),StrucParam.N_Tr-length(n2_set));
    for m=1:length(m_set)
        for r=1:1:StrucParam.N_Tr-length(n2_set)
            
            G_mr_s = zeros(length(m_set),1); % preallocation
            
            for s=1:length(m_set)
                G_mr_s(s) = (a_matr(m,s) * alpha_m(s) - a_sq(m,s) * eig2_m(length(n2_set) + r)) * F_mr_m(s,r);
            end
            
            G_mr_m(m,r) = (StrucParam.Z0/(2*pi/StrucParam.lambda)/StrucParam.eps2) * sum(G_mr_s);
            
        end
    end
    
else G_mr_m = [];
end

end