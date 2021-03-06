function [R_tot,T_tot,R_orders,T_orders,R,T,m_set,n1_set_ind,n2_set_ind,beta1_m_rot,beta2_m_rot,beta1_m_anti_rot,beta2_m_anti_rot,k0,eig1_p,  eig1_m,eig2_p,  eig2_m]=C_method_compute(StrucParam) %
%[m_set, n_set, alpha_m, beta1_m, beta2_m, n1_set, n1_set_ind, n2_set, n2_set_ind,alpha0,beta0] = CalculateAlphaBetaAndMNSets(StrucParam); % basic parameters
[m_set,  alpha_m, beta1_m, beta2_m,alpha0,beta1_ma, beta2_ma,beta0_a,beta1_m_rot,beta2_m_rot,beta1_m_anti_rot,beta2_m_anti_rot, beta0_rot_anti,n1_set, n1_set_ind, n2_set, n2_set_ind,k0]=CalculateAlphaBetaAndMNSets(StrucParam);
if strcmp(StrucParam.Fourier,'Numerical_int') || strcmp(StrucParam.Fourier,'Closed_form')     
   if strcmp(StrucParam.Fourier,'Closed_form') && isempty(StrucParam.Profile)
       disp('Warning: no closed_form expression for this profile, switch to numerical integration method instead');
       StrucParam.Fourier='Numerical_int';
   end
%% By Numerical integration or closed-form integration

%% Compute Fourier coefficient of a_diff
a_diff = ComputeDiff_of_a_array(StrucParam);

%% Evaluating L
if strcmp(StrucParam.Fourier,'Closed_form') && strcmp(StrucParam.Profile,'Harmonics')
    N_hm=length(StrucParam.An);
    if N_hm>1
        disp('Warning:there is no simple closed-form for L,switch to numerical integration instead');
    end
end

[L_mn_beta1, L_mk_beta2, L_mo_beta1] = Compute_L_arrays(m_set, StrucParam, beta1_m_rot,beta0_rot_anti, beta2_m_rot, n1_set, n1_set_ind, n2_set, n2_set_ind); %computation of L_mn, L_mk, and L_mo matrices

end
                                                                 
if strcmp(StrucParam.Fourier,'FFT')
%% By FFT
%% 2.Evaluate the FFT of function a_diff

a_diff_fun=@(x) eval(StrucParam.diff_a_x);
fm=StrucParam.fm;
a_diff_vec=F_series_gen(a_diff_fun,fm,StrucParam.dx2*StrucParam.Scale,StrucParam.N_Tr,StrucParam);

a_diff=[flipud(a_diff_vec);conj(a_diff_vec(2:end))];

%% Evaluating L
b_fun=@(x) eval(StrucParam.b_x);
[L_mo_beta1,L_mk_beta2,L_mn_beta1]=GenerateFFieldsChand(b_fun,beta0_rot_anti,StrucParam.N_Tr,StrucParam.dx2*StrucParam.Scale,m_set(1),m_set(end),n1_set ,n2_set , beta1_m_rot,beta2_m_rot,StrucParam);
end

 %% Eigvalue and eigvectors
[eig1, eig2, vect1, vect2] = SolveEigenvalues(a_diff, alpha_m, beta1_m, beta2_m, StrucParam,k0); %solve the eigenvalue problems in two media

%length(n1_set)
[eig1_p, vect1_p, eig1_m, ~] = SortEigenvaluesAndVectors(eig1, vect1,beta1_m_rot(n1_set_ind)/k0,beta1_m_anti_rot(n1_set_ind)/k0,StrucParam.kVecImagMin);
eig1_p=eig1_p*k0;
eig1_m=eig1_m*k0;

%length(n2_set)
[eig2_p, ~, eig2_m, vect2_m] = SortEigenvaluesAndVectors(eig2, vect2, beta2_m_anti_rot(n2_set_ind)/k0,beta2_m_rot(n2_set_ind)/k0,StrucParam.kVecImagMin);
eig2_m=eig2_m*k0;
eig2_p=eig2_p*k0;

%%{
if (StrucParam.measurement)
    clear eig1 eig2 vect1 vect2 beta1_m  beta2_m beta1_m_anti_rot beta2_m_anti_rot
end
%}

%% Assembling F
%% F_mn_R_p, F_mk_R_m, F_mo_R_in, F_mq_p, F_mr_m
F_mn_R_p = L_mn_beta1;
F_mk_R_m = L_mk_beta2;
F_mo_R_in = L_mo_beta1;
%F_mq_p = vect1_p(1:StrucParam.N_Tr,length(n1_set)+1:end);
F_mq_p = vect1_p(1:StrucParam.N_Tr,length(n1_set)+1:StrucParam.N_Tr);
%F_mr_m = vect2_m(1:StrucParam.N_Tr,length(n2_set)+1:end);
F_mr_m = vect2_m(1:StrucParam.N_Tr,length(n2_set)+1:StrucParam.N_Tr);

%% cut small elements, if the flag is set
if StrucParam.cut == 1
    F_mn_R_p = CutSmallArrayElements(F_mn_R_p, StrucParam.accuracy);
    F_mk_R_m = CutSmallArrayElements(F_mk_R_m, StrucParam.accuracy);
    F_mo_R_in = CutSmallArrayElements(F_mo_R_in, StrucParam.accuracy);
    F_mq_p = CutSmallArrayElements(F_mq_p, StrucParam.accuracy);
    F_mr_m = CutSmallArrayElements(F_mr_m, StrucParam.accuracy);
end

%% Assembling G
[G_mn_R_p, G_mk_R_m, G_mo_R_in, G_mq_p, G_mr_m] = Compute_G_matrices(m_set, n1_set, n1_set_ind, n2_set, n2_set_ind, ...
                                                                     alpha_m, a_diff, beta1_m_rot,beta0_rot_anti, beta2_m_rot,L_mn_beta1, L_mk_beta2, ...
                                                                     L_mo_beta1, eig1_p, F_mq_p, eig2_m, F_mr_m, StrucParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matching boundary conditions
GF_matrix = [F_mn_R_p F_mq_p -F_mk_R_m -F_mr_m;
                   G_mn_R_p G_mq_p -G_mk_R_m -G_mr_m];

GF_col = -[F_mo_R_in; G_mo_R_in];

if ~(StrucParam.n2==1i*inf)
Amplitudes = GF_matrix\GF_col;
else
    GF_matrix = [F_mn_R_p F_mq_p ;
                       G_mn_R_p G_mq_p];

   GF_col = -[F_mo_R_in; G_mo_R_in];
   
   Amplitude_plus=GF_matrix\GF_col;
   Amplitudes =[Amplitude_plus;Amplitude_plus*0];    
end

if strcmp(StrucParam.Pol,'BOTH')
    Polarization='TE';
    StrucParam = SetConstantsByPolarization(StrucParam,StrucParam.n1, StrucParam.n2, Polarization);
    [G_mn_R_p, G_mk_R_m, G_mo_R_in, G_mq_p, G_mr_m] = Compute_G_matrices(m_set, n1_set, n1_set_ind, n2_set, n2_set_ind, ...
                                                                     alpha_m, a_diff, beta1_m_rot,beta0_rot_anti, beta2_m_rot,L_mn_beta1, L_mk_beta2, ...
                                                                     L_mo_beta1, eig1_p, F_mq_p, eig2_m, F_mr_m, StrucParam);
     
     GF_matrix = [F_mn_R_p F_mq_p -F_mk_R_m -F_mr_m;
             G_mn_R_p G_mq_p -G_mk_R_m -G_mr_m];

    GF_col = -[F_mo_R_in; G_mo_R_in];
    
    
    if ~(StrucParam.n2==1i*inf)
    Amplitudes_TE = GF_matrix\GF_col; 
else
    GF_matrix = [F_mn_R_p F_mq_p ;
                       G_mn_R_p G_mq_p];

   GF_col = -[F_mo_R_in; G_mo_R_in];
   
   Amplitude_plus=GF_matrix\GF_col;
    Amplitudes_TE  =[Amplitude_plus;Amplitude_plus*0];    
end
    
    
    Polarization='BOTH';
    StrucParam = SetConstantsByPolarization(StrucParam,StrucParam.n1, StrucParam.n2, Polarization);  %Recover the polarization, for the next loop step.
end

%% calculation of reflection efficiencies
if strcmp(StrucParam.Pol,'BOTH')
   %fprintf('Reflection efficiency for %s polarization\n', StrucParam.Pol);   
   R_TM= zeros(StrucParam.N_Tr,1);
   R_TM(n1_set_ind)=real(beta1_ma(n1_set_ind)')./(-beta0_a).*abs(Amplitudes(1:length(n1_set)).^2);
   R_TE= zeros(StrucParam.N_Tr,1);
   R_TE(n1_set_ind)=real(beta1_ma(n1_set_ind)')./(-beta0_a).*abs(Amplitudes_TE(1:length(n1_set)).^2);
   R=[R_TM,R_TE];
   %fprintf(' First colume: orders; Second colume: for TM; Third colume: for TE\n');
   R_orders=[n1_set',R_TM(n1_set_ind),R_TE(n1_set_ind)];
   R_tot=[sum(R_TM),sum(R_TE)];
else
    %fprintf('Reflection efficiency for %s polarization\n', StrucParam.Pol);
R = zeros(StrucParam.N_Tr,1);
R(n1_set_ind)=real(beta1_ma(n1_set_ind)')./(-beta0_a).*abs(Amplitudes(1:length(n1_set)).^2);
R_orders=[n1_set',R(n1_set_ind)];
R_tot=sum(R);
end


%% calculation of transmission efficiencies
if strcmp(StrucParam.Pol,'BOTH')
   %fprintf('Transmission efficiency for %s polarization\n', StrucParam.Pol);
   T_TM= zeros(StrucParam.N_Tr,1);
   T_TE= zeros(StrucParam.N_Tr,1);
   if ~isempty(n2_set)
   T_TM(n2_set_ind)=StrucParam.eps1*(-1)*beta2_ma(n2_set_ind)'./(StrucParam.eps2*(-beta0_a)).*abs(Amplitudes(StrucParam.N_Tr+1:StrucParam.N_Tr+length(n2_set_ind))).^2;
    Polarization='TE';
    StrucParam = SetConstantsByPolarization(StrucParam,StrucParam.n1, StrucParam.n2, Polarization);
   T_TE(n2_set_ind)=StrucParam.eps1*(-1)*beta2_ma(n2_set_ind)'./(StrucParam.eps2*(-beta0_a)).*abs(Amplitudes_TE(StrucParam.N_Tr+1:StrucParam.N_Tr+length(n2_set_ind))).^2;
   Polarization='BOTH';
   StrucParam = SetConstantsByPolarization(StrucParam,StrucParam.n1, StrucParam.n2, Polarization); 
   end
   T=[T_TM,T_TE];
   %fprintf(' First colume: orders; Second colume: for TM; Third colume: for TE\n');
   T_orders=[n2_set',T_TM(n2_set_ind),T_TE(n2_set_ind)];
   T_tot=[sum(T_TM),sum(T_TE)];
else
     %fprintf('Transmission efficiency for %s polarization\n', StrucParam.Pol);
T = zeros(StrucParam.N_Tr,1);
if ~isempty(n2_set)
T(n2_set_ind)= StrucParam.eps1*(-1)*beta2_ma(n2_set_ind)'./(StrucParam.eps2*(-beta0_a)).*abs(Amplitudes(StrucParam.N_Tr+1:StrucParam.N_Tr+length(n2_set_ind))).^2;
end
T_tot=sum(T);
T_orders=[n2_set',T(n2_set_ind)];
end

%% Plot field distribution if required.
%It is suggested to be only for single point, at most one wavelength.
%%{
     if strcmp(StrucParam.plot_fields,'YES') && ~(StrucParam.measurement)
        [imag_eig1p,imag_eig2n,imag_Vec1p,imag_Vec2n]=SortEigenvaluesChand(eig1,vect1,eig2,vect2,StrucParam.kVecImagMin,StrucParam.N_Tr);
        imag_eig1p=imag_eig1p*k0;
        imag_eig2n=imag_eig2n*k0;
        b_fun=@(x) eval(StrucParam.b_x); 
        if strcmp(StrucParam.Pol,'BOTH')
        StrucParam.Pol='TM';
        Plot_Intensity(alpha0,beta0_rot_anti,b_fun,Amplitudes,n1_set, beta1_m_rot, alpha_m, n2_set,imag_Vec1p, imag_eig1p, beta2_m_rot, imag_Vec2n,imag_eig2n, m_set(1), m_set(end), StrucParam);
        StrucParam.Pol='TE';
        Plot_Intensity(alpha0,beta0_rot_anti,b_fun,Amplitudes_TE,n1_set, beta1_m_rot, alpha_m, n2_set,imag_Vec1p, imag_eig1p, beta2_m_rot, imag_Vec2n,imag_eig2n, m_set(1), m_set(end), StrucParam);
       StrucParam.Pol='BOTH';
        else
         Plot_Intensity(alpha0,beta0_rot_anti,b_fun,Amplitudes,n1_set, beta1_m_rot, alpha_m, n2_set,imag_Vec1p, imag_eig1p, beta2_m_rot, imag_Vec2n,imag_eig2n, m_set(1), m_set(end), StrucParam);
        end
    end
      %}

end