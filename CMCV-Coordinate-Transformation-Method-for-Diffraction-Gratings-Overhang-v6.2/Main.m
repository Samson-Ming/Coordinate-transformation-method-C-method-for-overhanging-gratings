clear; 
clc ; 
close all;

%% Seting initial parameters in 'SetInitialParameters.m' first, and then run 'Main.m'.
%Input the initial parameters in the SetInitialParameters.m, and import them here.
SetInitialParameters

%% If you still want to change some parameters, you can set them here or later, mainly by updating the struct StrucParam.
%% BE CAREFUL with dispersion and polarization, 
%% if changed, then update with StrucParam = SetConstantsByPolarization(StrucParam,n1, n2, Polarization);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start calculations
%% Set the measurement needed below:
%-------------------------------------------------------     
% setup measurement
% template measurement
% 0 --- shows all diffraction efficiencies and total reflection and transmission
% 1 --- dependence of the diffraction efficiencies on the trunction number
% 2 --- dependence of the diffraction efficiencies on the wavelength. Update StrucParam if dispersive.
% 3 --- dependence of the diffraction efficiencies on the incident angle theta
% 4 --- dependence of the diffraction efficiencies on the refractive index $n_1$ ($n_g$). Update StrucParam if dispersive.
% 5 --- dependence of the diffraction efficiencies on the refractive index $n_2$ ($n_r$). Update StrucParam if dispersive.
% 6 --- dependence of the diffraction efficiencies on the grating period dx
% 7 --- dependence of the diffraction efficiencies on the rotation angle Phi
% 8 --- dependence of the diffraction efficiencies on the user defined parameters, 
%        especially for optimization design. Define any of the other measurement here.

StrucParam.measurement=0;

% !!! setup of measured diff. order when measurement>0
diffraction_efficiencies_c=1; %1--R, 2--T.M
Pol_ind=1; %1--TM,2--TE; Used when StrucParam.Pol='BOTH' in R and T;
if ~strcmp(StrucParam.Pol,'BOTH')
    Pol_ind=1;
end
studying_order=-1;  %studying_order must in m_set. Otherwise set to 1;

%% Scale factor, will change when sweeping wavelength or sweeping period
%% Several choice of StrucParam.Scale
%StrucParam.Scale=1;  %No Nomalization
%StrucParam.Scale=1/StrucParam.dx; %Nomalize the period to 1 (no unit).
%StrucParam.Scale=2*pi/StrucParam.wavelength;  %Nomalize the wavenumber to 1 (no unit)
StrucParam.Scale=1/StrucParam.wavelength; %Nomalize the wavelength to 1 (no unit).
StrucParam.lambda =StrucParam.wavelength*StrucParam.Scale; %wavelength need scaling here.

tic;
switch StrucParam.measurement
    case 0
       %% Single point Diffraction efficiency              
    case 1
    %% Convergence test, sweep N_Tr
    N=110:1:140;
    
    case 2
      %% Update StrucParam.Scale if wavelength or wavenumber is normalized.
      %% Update StrucParam if dispersive.
        Wavelength=linspace(0.3,0.8,20);
        
        %{
        %Dispersion: Sellmier
        A=1;
        B=[1.4313493 0.65054713 5.3414021];
        C=[0.0726631 0.1193242 18.028251];
        n2=sellmeier(Wavelength,A,B,C);
        %}
        
        %{
        %Dispersion: Drude
        omega=2*pi*nm2Hz(Wavelength*1e3);  %omega in [rad/s], usually wavelength is in [um]
        n2=drude(omega,1,1.32e16,1.2e14);
        %}
        
        %{
        %Dispersion: interpolation
        n2=rix_spline(Wavelength,'gold_palik.txt');
        %}

                   
    case 3
        % 3 --- dependence of the diffraction efficiencies on the incident angle theta
        Theta=linspace(0,50,51)*pi/180;  %unit: [rad].(
        
    case 4
        %% Update StrucParam.
         n1=linspace(1,1.5,11);
        
    case 5
        %% Update StrucParam.
         n2=linspace(1.5,2.5,11);   %for complex n2, define the real part and complex part seperately.

    case 6
        %% Inapplicable for grating=0,5,6,and maybe 7.
        %% Update StrucParam.Scale if period is normalized.
          Period=linspace(0.4,0.65,20);
        
    case 7
        %% ONLY some profile disribed in the ORIGINAL coordinate (grating=0,2,3,5,6, and maybe 7)
        %% There might be problem if the a(x) is not function after rotation.
        %% If Phi is in the correct range, this curve should be a flat line 
        Phi=linspace(28,32,5)*pi/180;  %unit: [rad].
       
    case 8
        %% DO NOT use "index" as iteration indicator because it has alread been used in SetInitialParameters.m
        %% If you change the grating profile parameters, remember to run "setup_grating_profile;" to update the grating profile.
        %% Put the sweep script in a M file. Here is a typlical example below.
        
         measurement_file='sweep_grating_height.m';
         
         %% Also you can define you own main file elsewhere, and do not use the default "Measurement" here, 
               %just use the "SetInitialParameters" and function "C_method_compute" as a blackbox    
end
Measurements

%save metal_deep_overhanging_power_sine_convergence_PPT.mat N Convergence R_convergence StrucParam

%{
save metal_deep_overhanging_power_sine_convergence_PPT_N_Tr=401_TM.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal
%}

%{
save Plumey1997_LD_right_triangle_NO_PPT_N_Tr=27_1percent_TE.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal
%}

%{
R_ref=[0.266477888803202,0.247755429397570,0.141235922156067];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Metal_right_triangle_NO_PPT_convergence_Phi0_TE2.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy
   
%}

%{
R_ref=[0.464534519597190   0.106485838427399   0.015466816023207];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save LD_Metal_right_triangle_NO_PPT_convergence_Phi0_TM.mat   
%}

%{
save Dielectric_overhang_trapezoid_NO_PPT_N_Tr=299_1percent_TE.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal
%}

%{
R_ref=[0.0116110339181890	0.00397254128892009	0.00269456222006578];
T_ref=[0.0984302860119821	0.476629669959978	0.0994663070054246	0.268814629169694	0.0384145463940903];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
%save Dielectric_overhang_trapezoid_NO_PPT_convergence_Phi0_TE.mat N Convergence R_convergence T_convergence StrucParam...
%        RT_convergence RT_ref Accuracy
   
%}

%{
save Dielectric_overhang_trapezoid_NO_PPT_N_Tr=299_1percent_TM.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal
%}

%{
R_ref=[];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Dielectric_overhang_trapezoid_NO_PPT_convergence_Phi0_TM.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy
   
%}

%{
save Metal_overhang_sinusodial_h=0.5lam_Phi=60_PPT=7_N_Tr=13_1percent_TE.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal
%}

%{
R_ref=[0.344941437878085,0.534255315313941];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Metal_overhang_sinusodial_h=0.5lam_Phi=60_PPT=7_convergence_Phi0_TE.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy
   
%}

%{
save Metal_overhang_sinusodial_PPT=7_h=0.5lam_Phi=60_N_Tr=19_1percent_TM.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal                                                                                                      
%}
%save Metal_overhang_sinusodial_PPT=7_h=0.5lam_Phi=60_N_Tr=19_1percent_TM_Near_filed.mat     
%{
R_ref=[0.753199898155597,0.0400481471266204];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Metal_overhang_sinusodial_PPT_PPT=7_h=0.5lam_Phi=60_convergence_TM.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy   
%}

%{
save Metal_overhang_sinusodial_PPT=7_h=lam_Phi=60_N_Tr=21_1percent_TM.mat R_tot T_tot R_orders T_orders R T ...
                                                                                                                m_set n1_set_ind n2_set_ind...
                                                                                                                beta1_m_rot beta2_m_rot beta1_m_anti_rot beta2_m_anti_rot k0...
                                                                                                                eig1_p eig1_m eig2_p eig2_m...
                                                                                                                timerVal                                                                                                      
%}

%{
R_ref=[0.726026207305534,0.0820952701445662];
T_ref=[];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Metal_overhang_sinusodial_PPT_PPT=7_h=lam_Phi=60_convergence_TM.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy   
%}

%{
save Dielectric_overhang_trapezoid_half_h_NO_PPT_N_Tr=105_1percent_TE.mat
%}

%{
R_ref=[0.00795593990107011,0.00174225525041201,0.00791078068522617];
T_ref=[0.0425553602381309,0.334609217015269,0.310258547411985,0.281629835244023,0.0133448942078767];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Dielectric_overhang_trapezoid_half_h_NO_PPT_convergence_Phi0_TE.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy StrucParam
   
%}

%{
save Dielectric_overhang_trapezoid_half_h_NO_PPT_N_Tr=257_1percent_TM.mat 
%}

%{
R_ref=[0.0107100247626546,0.000404796285108852,0.00440433589164716];
T_ref=[0.0115596868078754,0.322548395911140,0.465857399052242,0.178003022159935,0.00652128980189012];
RT_ref=[R_ref,T_ref]';
RT_convergence=[R_convergence;T_convergence];
%Accuracy=abs(RT_convergence-RT_ref)./RT_ref;
Accuracy=max(abs(RT_convergence-RT_ref)./RT_ref);
save Dielectric_overhang_trapezoid_half_h_NO_PPT_convergence_Phi0_TM.mat N Convergence R_convergence T_convergence StrucParam...
        RT_convergence RT_ref Accuracy StrucParam
   
%}

%save metal_deep_overhanging_power_sine_convergence_PPT=9_k-3_MP_N=1001_BOTH.mat
%save metal_deep_overhanging_power_sine_convergence_PPT=9_LP_convergence_k-3_TM.mat