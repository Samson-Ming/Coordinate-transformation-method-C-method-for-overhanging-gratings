function [L_mn_beta1,L_mk_beta2,L_mo_beta1]=Compute_L_arrays(m_set,StrucParam,beta1_m,beta0,beta2_m,n1_set,n1_set_ind,n2_set,n2_set_ind)
idx_m = 1;
idx_n = 1;
%% L_mn_beta1
if ~isempty(n1_set)
    L_mn_beta1 = zeros(length(m_set), length(n1_set));% preallocation
    for m = m_set(1):1:m_set(end)
        for n = n1_set(1):1:n1_set(end)
            L_mn_beta1(idx_m, idx_n) = L_eval(beta1_m(n1_set_ind(idx_n)), m-n, StrucParam);
            idx_n = idx_n+1;
        end
        idx_n = 1;
        idx_m = idx_m + 1;
    end
else L_mn_beta1 = [];
end

%% L_mk_beta2
idx_m = 1;
idx_k = 1;
if ~isempty(n2_set)
    L_mk_beta2 = zeros(length(m_set), length(n2_set));
    for m = m_set(1):1:m_set(end),
        for k = n2_set(1):1:n2_set(end),
            L_mk_beta2(idx_m, idx_k) = L_eval(beta2_m(n2_set_ind(idx_k)), m-k, StrucParam);
            idx_k = idx_k + 1;
        end
        idx_k = 1;
        idx_m = idx_m + 1;
    end
else L_mk_beta2 = [];
end

%% L_mo_beta0
idx = 1;
L_mo_beta1 = zeros(length(m_set),1); %preallocation
for m = m_set(1):1:m_set(end)
    L_mo_beta1(idx, 1) = L_eval(beta0, m, StrucParam);
    idx = idx + 1;
end

end


