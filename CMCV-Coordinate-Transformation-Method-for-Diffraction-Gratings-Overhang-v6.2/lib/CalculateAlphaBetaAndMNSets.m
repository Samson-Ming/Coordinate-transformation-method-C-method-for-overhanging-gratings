function [m_set,alpha_m, beta1_m, beta2_m,alpha0,beta1_ma, beta2_ma,beta0_a,beta1_m_rot,beta2_m_rot,beta1_m_anti_rot,beta2_m_anti_rot, beta0_rot_anti,n1_set, n1_set_ind, n2_set, n2_set_ind,k0]=CalculateAlphaBetaAndMNSets(StrucParam)
k0=2*pi/StrucParam.lambda;

%% ----------------------------In the ORIGINAL coodinate----------------------------------------------------
alpha0_a=StrucParam.n1*k0*sin(StrucParam.theta);
K_a = 2*pi/(StrucParam.dx*StrucParam.Scale); %inverse lattice vector
m1 = -floor(alpha0_a/K_a)-((StrucParam.N_Tr-1)/2);
m2 =-floor(alpha0_a/K_a)+((StrucParam.N_Tr-1)/2);
StrucParam.N_Tr=m2-m1+1;

m_set = m1:1:m2; % array of m according to Li
n_set = 1-StrucParam.N_Tr:1:StrucParam.N_Tr-1; % array of n according to Li

%% ----------------------------In the ORIGINAL coodinate----------------------------------------------------
alpha_ma = alpha0_a+K_a*m_set;
beta1_ma = sqrt((StrucParam.n1*k0)^2-alpha_ma.^2);
beta0_a=-beta1_ma(1-m1);
beta2_ma = -sqrt((StrucParam.n2*k0)^2-alpha_ma.^2);

%%  ----------------------------In the ROTATED coodinate----------------------------------------------------
Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];
%inv_Rot=[cos(StrucParam.Phi),sin(StrucParam.Phi);-sin(StrucParam.Phi),cos(StrucParam.Phi)];

%[alpha_m;beta2_m]
alpha_beta2_rot=Rot*[alpha_ma;beta2_ma];
beta2_m_rot=alpha_beta2_rot(2,:);

alpha_beta2_anti_rot=Rot*[alpha_ma;-beta2_ma];
beta2_m_anti_rot=alpha_beta2_anti_rot(2,:);

%[alpha_m;beta1_m]
alpha_beta1_rot=Rot*[alpha_ma;beta1_ma];
beta1_m_rot=alpha_beta1_rot(2,:);

alpha_beta1_anti_rot=Rot*[alpha_ma;-beta1_ma];
beta1_m_anti_rot=alpha_beta1_anti_rot(2,:);
beta0_rot_anti=beta1_m_anti_rot(1-m1);

%%  ----------------------------For the New grating in the ROTATED coodinate----------------------------------------------------
%{
thetap=asin(sin(StrucParam.theta)/cos(StrucParam.Phi));
thetaa = StrucParam.theta+StrucParam.Phi;
theta_rot = StrucParam.theta+StrucParam.Phi;
thetap=asin(sin(StrucParam.theta)/cos(StrucParam.Phi));
%}

alpha_m=alpha_ma/cos(StrucParam.Phi);
alpha0=alpha_m(1-m1);
beta1_m = sqrt((StrucParam.n1*k0)^2-alpha_m.^2);
beta0=-beta1_m(1-m1);
beta2_m = -sqrt((StrucParam.n2*k0)^2-alpha_m.^2);

%When some element of beta1_m or beta2_m is equal to zero, M are singular,add a small increment to continue
%When beta1_ma or beta2_ma is singular, i.e. at Rayleigh anomoly, the symstem may also be singular
if (min(abs(beta1_m)<StrucParam.tol) || min(abs(beta2_m)<StrucParam.tol)) ||  (min(abs(beta1_ma)<StrucParam.tol) || min(abs(beta2_ma)<StrucParam.tol))              
    disp('Singularity, change the incident angle slightly');
    StrucParam.theta =StrucParam.theta + 1*pi/180000;
    
    alpha0_a=StrucParam.n1*k0*sin(StrucParam.theta);
    alpha_ma = alpha0_a+K_a*m_set;
    beta1_ma = sqrt((StrucParam.n1*k0)^2-alpha_ma.^2);
    beta0_a=-beta1_ma(1-m1);
    beta2_ma = -sqrt((StrucParam.n2*k0)^2-alpha_ma.^2);
    
    %[alpha_m;beta2_m]
    alpha_beta2_rot=Rot*[alpha_ma;beta2_ma];
    beta2_m_rot=alpha_beta2_rot(2,:);
    alpha_beta2_anti_rot=Rot*[alpha_ma;-beta2_ma];
    beta2_m_anti_rot=alpha_beta2_anti_rot(2,:);   
    
    %[alpha_m;beta1_m]
    alpha_beta1_rot=Rot*[alpha_ma;beta1_ma];
    beta1_m_rot=alpha_beta1_rot(2,:);
    alpha_beta1_anti_rot=Rot*[alpha_ma;-beta1_ma];
    beta1_m_anti_rot=alpha_beta1_anti_rot(2,:); 
    beta0_rot_anti=beta1_m_anti_rot(1-m1);
        
   %%-------------------------------------
   alpha_m=alpha_ma/cos(StrucParam.Phi);
   alpha0=alpha_m(1-m1);
   beta1_m = sqrt((StrucParam.n1*k0)^2-alpha_m.^2);
   beta0=-beta1_m(1-m1);
   beta2_m = -sqrt((StrucParam.n2*k0)^2-alpha_m.^2);
end




%Sort out Propagating modes, either use beta_m_rot or beta_ma, they are equivalent.
idx=find(~imag(beta1_ma) & beta1_ma~=0);
n1_set = m_set(idx);
n1_set_ind = idx.';
idx=find(~imag(beta2_ma) & beta2_ma~=0);
n2_set = m_set(idx);
n2_set_ind = idx.';

end