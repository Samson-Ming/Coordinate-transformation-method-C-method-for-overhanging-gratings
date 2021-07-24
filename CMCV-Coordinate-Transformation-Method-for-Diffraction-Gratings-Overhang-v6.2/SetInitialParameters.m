%% Seting initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All of the intial parameters are stored in the structure "StrucParam". %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% You can visit most parameters by accessing StrucParam  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% and these these parameters by updating StrucParam  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% First construct an empty StrucParam
format long;
addpath('measurements');
addpath('grating profiles');
addpath('dispersion models');
addpath('lib');
StrucParam=struct();  %Plesea do not change this row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Control parameters:
%% Check intial parameters on not?
StrucParam.Check_parameters='NO';    %'NO', default, not check initial parameters.
                                     %'YES', check and display intial parameters, recommended at the begining;                                                          
 %% flag - accuracy rounding 
StrucParam.cut = 0;                 % 0 - default, no rounding, 1 - round to zero small elements in the array a_diff and all F and G matrices
StrucParam.accuracy = 1e-15;   % set the threshold for rounding

%% accuracy of imaginary part of k-vector
StrucParam.kVecImagMin = 1e-12;% default, a k-vector is set to be real if abs(Imag(k)) is smaller than this threshold
StrucParam.tol = 1e-15;% default, abs(beta1_m or beta2_m) smaller than this threshold is regarded to be close to Rayleigh anomoly

%% truncation order of the harmonics
%StrucParam.N_Tr = 2*15 + 1; % odd number
N=200;
StrucParam.N_Tr = 2*N + 1; % This is just the initial set of trunction number, the actual trunction number will be adjusted in the program.
                                          % type StrucParam.N_Tr to check the real trunction number.
                                          
%% Very important: Data Precision transform for numerical stability at LARGE StrucParam.N_Tr.
%% In fact it is very difficult to decide whether to use or not the Data Precision Transform technique.
%% The best way to determine StrucParam.Precision_transform is by converge test.
%% and by set to 0 if it reports ERROR or if the ultimate results are apprarently WRONG!
%% Also the results need verifying with converged result by OTHER method (FMM/DM/LPEM etc.)
StrucParam.PPT='YES';  %'NO':default, not apply the perturbative preconditioning technique.
                                  %'YES': apply the perturbative preconditioning technique.
                                  % Specially for smooth profile,like sinusodial grating.
                                  % Highly recommended for deep smooth lossy grating with large trunction number.
                                  % For profile with sharp corners,like triangle grating, may cause strange unstabitlity, then change to 'NO'.
StrucParam.PPT_method='LP';
                           %Default,'LP',low precision, chop function;
                           %'MP',multiple precision, Advanpix. 
                           % 'MP' is ONLY valid in Windows AND you MUST have purchased the Advanpix toolbox, and copy it to the 'Advanpix' filefolder.
                           %'Sym', Symbolic Math toolbox;
                           %'Round', use the simple round function;
if strcmp(StrucParam.PPT_method,'MP' ) 
    addpath('Advanpix');
end
                                  
StrucParam.decimal_digits=7;  % Default. Other recommended:decimal_digits=8,9.
if strcmp(StrucParam.PPT,'YES')
StrucParam.kVecImagMin = 10^-(StrucParam.decimal_digits-3);%  a k-vector is set to be real if abs(Imag(k)) is smaller than this threshold
%StrucParam.kVecImagMin = 1e-12; % default,
                                 % StrucParam.kVecImagMin is quite intricate, need to adjust for different problem and StrucParam.decimal_digits
end
%% flag - Method used to calculate the Fourier coefficients.
StrucParam.Fourier='FFT';  
                                      %'FFT': Default,use FFT to calculate the Fourier coefficients.
                                      %'Closed_form': use closed-form expression. Always use closed form expression if applicable.
                                      %'Numerical_int':  use numerical integration algorthm (like quadgk) to calculate Fourier coefficients, usually avoid this method.                                                                                           
%% FFT orders, used when StrucParam.Fourier='FFT';
StrucParam.fm = 14;   %Numbers of sampling point 2^fm+1, i.e. FFT orders: fm;
% Detailed comparison of differnt methods to calculate the Fourier coefficients.
%%%%%%%%%%%%%%%%%%%%%%%% Compare between the three methods to calculate Fourier coefficients  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                     
%{
                              Generity                                  Speed                     Stability                          Precision

FFT                          General                                  medium                   stable                         relatively low

Closed_form          Several profiles                            fast                       stable                              high

Numerical_int       Best for explicit profiles                 slow              sometimes unstable             relative high

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                               
%% Incident light
%% angle of incidence -> [0, pi/2), not including pi/2 - grazing angle
%%  --------------------------------------------------------------------------------
% BE CAREFUL: This is the incident angle in the original UN-Rotated coordinate
%StrucParam.theta = 45*pi/180; 
%StrucParam.theta = 25*pi/180; 
StrucParam.theta =asin(1/3); 
%StrucParam.theta = 0*pi/180;   %normal incidence. 
                                                 %For theta=0, sometimes
                                                 %errors can occur due tosingularity of the eigenMatrix. Littrow mounting?
                                                 %If such situation occurs, just add a tiny value to theta to walkaround the numerical difficulty, 
                                                 %say set StrucParam.theta = 0*pi/180+1e-3*pi/180;
%% wavelength in any unit, recommend in micron for visible range.
StrucParam.wavelength=0.6328;   %unit: a.u.

%% polarization: 'TE' , 'TM' or 'BOTH'
Polarization='TM'; %Current only for collinear mounting, not for conical mounting. For 'BOTH', the calculation time only increases a little.
   
%% Grating properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scale factor: important.
%Scale all of the geometrical parameters (length quantities) to a suitable scale, 
%in order to gurantee the numerical precision, especially in the integration
% Note that Maxwell's is scale invarant, thus the length quantity can be in arbitary unit (a.u.).
StrucParam.Scale = 1;  
                       %DO NOT multiple the StrucParam.Scale when inputing the structure parameters here.
                       %BUT Use StrucParam.Scale in the expression of a_x, diff_a_x, and other function part when length quantities are scaled, 
                       %like ComputeDiff_of_a_array.m,L_eval.m,GenerateFFieldsChand.m and CalculateAlphaBetaAndMNSets.m
                       %just add '*StrucParam.Scale' to each geometrical length parameters (only for length quantities),
                       %NOTE that NO scale is needed for all the angle related quantities (explicite phase given by a number,slope, incident angle, azimuth angle, polarization angle...).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% material properties
% refractive indicies (sign for metals n_r+1i*n_i)
n1 = 1;       % refractive index of the first material (ambient), also the incident media. Must be lossless (real number).

%n2 =sqrt(-21+60.4i);  % Metal, refractive index of the second material, can be lossy (complex number) or lossless (real or pure imaginary number).
n2 =0.2+6i;
%n2 =1+5i;
%n2 = 1.5;  % Dielectric;
%n2 = 1i*inf;  % Perfect conductor 

%% arrangement for TE and TM polarizations
StrucParam = SetConstantsByPolarization(StrucParam,n1, n2, Polarization); 

%%%%%%%%%%%%%%%%%%%  Each time if n1,n2 or Polarization have been changed, %%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%  StrucParam must be updated by this commond.          %%%%%%%%%%%%%%%%%%%

%% grating profile type;
StrucParam.Profile=[]; % Leave it empty here, and will define later in grating type.
                                % Only valid for several kinds of profiles, and leave empty for general type. 

                                                                  
%% Plot field distribution if required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%  It is suggested to be only for single point, at most one   wavelength because it is quite time consuming.%%%%%%%%%%%%%%%%%%%%%%%%%%% 
StrucParam.plot_fields='NO'; %'NO',default,not calculate the field. 'YES': calculate the field map; It is suggested to be only for single point (measurement=0), at most one wavelength.
StrucParam.save_plot='NO';  %'NO',default,not save the field map; 'YES': save the field map to a png file.
StrucParam.Num_of_Period=1; %Plot field in several numbers of periods;
StrucParam.upper=0.1; %h_upper=StrucParam.upper*h_grating;
StrucParam.lower=0.1; %h_lower=StrucParam.lower*h_grating;
StrucParam.point_x=200;  %Sampling point number in the periodic direction IN PER PERIOD; Large number can significantly increase computation time.
                                      %TOTAL Number in the period direction=StrucParam.Num_of_Period*StrucParam.point_x
StrucParam.point_y=200;  %Sampling point number in the direction perpendicular to grating. Large number can significantly increase computation time.

%% Define grating profile

%% strucutre geometric parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Note: for diff_a_x, at discontinuity, diff_a_x=(diff_a_x(left_limit)+diff_a_x(rigth_limit))/2;
%---------------------------------------------------
% Setup of the grating
%  grating profile=b(x)+tan(Phi)*x;
% b(x)=
%  0 --- General Polyline grating in the ORINGINAL coordinate, Type 1.
%  1 --- General Polyline grating in the ROTATED coordinate, Type 2.
%  2 --- Full triangle grating in the ORIGINAL coordiante.
%  3--- Full trapezoid grating in the ORIGINAL coordiante.
%  4 --- Harmonic grating in the ROTATED coordinate, type1, explicite. (described by sumation of forier harmonics)
%  5 --- Harmonic grating in the ORIGINAL coordiante, type2, implicite, fit dicreticized sampling points by sumation of forier harmonics. 
%  6 --- Pure sampling points in the ORIGINAL coordinate, not fit by any function, use FFT. Curently not very good.
%  7 --- Power-sine in the ROTATED coordinate.
%  8 --- User defined grating from file (in the grating profiles filefolder).

grating=7; % switch grating
StrucParam.dx=1; % grating period in the ORIGINAL coordiante [a.u.] (exception grating=5,6)

%% Below are mainly for a(x) if not specified elsewhere.
% the detailed definition can be found in the grating profiles filefolder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%  BE CAREFUL with the ORIGINAL coordiante or ROTATED coordinate.  %%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch grating
    case 0
        %StrucParam.Phi=acot(0.36/(0.4211-0.425/2)); %Rotation angle, [rad]. 
        %StrucParam.Phi=40*pi/180; %Rotation angle, [rad]. 
        StrucParam.Phi=pi/2-atan(2);
      %% In the ORIGINAL coordiante, define the general Polyline grating by the coordinate of each node.
        %The first one must (x,y)=(0,0), and the last one must be (x,y)=(StrucParam.dx,0);
        %MUST be in row vector, nodes should not overlap with each other.
        %StrucParam.xn_ORI=[0 0.0618 0.4211 0.3632 StrucParam.dx];  %This is an overhang partial triangle example.
        %StrucParam.yn_ORI=[0 0 0.36 0 0];
        
        StrucParam.xn_ORI=[0 1/4 5/8 7/8 3/4 1]*StrucParam.dx;
        StrucParam.yn_ORI=[0   0     1    1     0   0]*StrucParam.dx/2;
                       
    case 1
        StrucParam.Phi=30*pi/180; %Rotation angle, [rad].
        StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate.
        
       %% In the ROTATED coordiante, define the general Polyline grating by the coordinate of each vertex.
        %The first one must (x,y)=(0,0), and the last one must be (x,y)=(StrucParam.dx2,0);
        %xn(m)<xn(m+1);        
        StrucParam.xn=linspace(0,StrucParam.dx2,100);
        StrucParam.yn=(0.360/2/cos(30*pi/180))*(1-cos(1*2*pi/StrucParam.dx2.*StrucParam.xn));  %This is an approximation of a sinusodial function by polylines, and it WORKS.
            
    case 2
        StrucParam.Phi=acot(1/(1-1/2)); %Rotation angle, [rad].
        
        %% In the ORIGINAL coordiante, the coordinate of the top vertex (x,y)=(c,h). 
        %Note that the other two vertexes have already been assumed as (0,0) and (StrucParam.dx,0);       
        StrucParam.c_ORI=1;  %This is an right triangle example.
        StrucParam.h_ORI=1;       
             
    case 3
         StrucParam.Phi=35*pi/180; %Rotation angle, [rad].

       %% In the ORIGINAL coordiante, the coordinate of the top two vertex (x,y)=(c1,h) and (x,y)=(c2,h);
        %Note that the other two vertexes have already been assumed as (0,0) and (StrucParam.dx,0); 
        StrucParam.c1_ORI=0.3;  %This is an overhang partial triangle example.
        StrucParam.c2_ORI=0.6;  %This is an overhang partial triangle example.
        StrucParam.h_ORI=0.26; 
    
   case 4
       StrucParam.Phi=60*pi/180; %Rotation angle, [rad]. 
       
      %% In the ROTATED coordiante, a(x') is represented by summation of COSINE functions
       % a(x')=sum_n=1:Nh_(A_n*cos(2*n*pi/dx2*x'+phi_n))
       %NOTE the DC component (i.e., the constan A_0*cos(phi_0)) is negnected
       %because it is merely a translation in the y' direction, which has no effect on the efficiency
       
       %{
       StrucParam.A0=0; %Constant component, 0-th frequency.
       StrucParam.An=[0.1 0.02 0 0 ]*0.525; % Fourier-coefficients amplitude     
        StrucParam.phi_n=[0.00 -5/9 0 0]*pi; % Fourier-coefficients phase
       %}
       
       %{
        StrucParam.A0=0; %Constant component, 0-th frequency.
       StrucParam.An=0.360/2/cos(30*pi/180); % Fourier-coefficients amplitude     
       StrucParam.phi_n=0*pi; % Fourier-coefficients phase
       %}
       
       StrucParam.A0=1/2/cos(60*pi/180);%Constant component, 0-th frequency.
       StrucParam.An=-1/2/cos(60*pi/180); % Fourier-coefficients amplitude     
       StrucParam.phi_n=0*pi; % Fourier-coefficients phase  
            
   case 5
       StrucParam.Phi=30*pi/180; %Rotation angle, [rad].  
      
      %% In the ORIGINAL coordiante, discrete the grating profile.
        % Usually sample the grating profile and save it in a text file, and load it here.
        %The first one must (x,y)=(0,0), and the last one must be (x,y)=(StrucParam.dx,0);
        %MUST be in row vector, nodes should not overlap with each other.
        input_sampling_file='overhang sinusodial grating.txt';
        
        %Need resampling to keep it uniformly-spaced.
        StrucParam.resample_order=14;     %Usually for convinience of FFT, N_resample=2^resample_order;
        
        %Numbers of COSINE functions used to appoximate the resampled data.
        StrucParam.Nh=1;  %NOTE NOT including the DC component (0 frequency), only the numbers of higher harmonics beginning from 1st order. 
        
        %Method to fit the resampled nodes to Fourier harmonics in ROTATED coordinate
        StrucParam.fit_method='FFT';    %'FFT', use FFT to retrieve Fourier coefficients, more accurate.
                                                      %'Inverse', default,use pseudo-inverse to retrieve Fourier coefficients, faster.
                                                                                                                        
    case 6
       StrucParam.Phi=30*pi/180; %Rotation angle, [rad].  
       
      %% Better to set StrucParam.Precision_transform=0; StrucParam.Fourier='FFT';  StrucParam.fm = Same order with sampling; 
       
      %% In the ORIGINAL coordiante, discrete the grating profile.
        % Usually sample the grating profile and save it in a text file, and load it here.
        %The first one must (x,y)=(0,0), and the last one must be (x,y)=(StrucParam.dx,0);
        %MUST be in row vector, nodes should not overlap with each other.
        %Need resampling to keep it uniformly-spaced.
        input_sampling_file='overhang sinusodial grating.txt';
        
        %Need resampling to keep it uniformly-spaced.
        StrucParam.resample_order=12;     %Usually for convinience of FFT, N_resample=2^resample_order;
        
    case 7
        %a(x_rot)=(h/cos(Phi))*[1-(1-sin(pi/(d*cos(Phi))*x_rot)^2alpha)^beta];
        StrucParam.Phi=30*pi/180; %Rotation angle, [rad].  
        StrucParam.h=10;   %Height of the grating in the ORIGINAL coordinate;
        %StrucParam.h=1;   %Height of the grating in the ORIGINAL coordinate;        
        StrucParam.Alpha=3;
        StrucParam.Beta=4;
        
        %{
        N=401;
                                    Rtot                              R-2                                    Time(s)
        LP             0.679114088634636    0.031599678167047            1.568992e+01
        MP            0.679106581513191    0.031600646869347            1.498928e+1(About)
        Sym          0.679126607582405    0.031599831369391             2.574403e+02
        Round       0.679095544459448    0.031608591416988             1.503884e+01
        %}
        
        %N=401,R_tot=0.680322894142198;
        %N=601,R_tot=0.679561545473441;
        %N=801,R_tot=0.679233480837092;
        %N=1001,R_tot=0.679406041486129,R-2=0.03159
       
    case 8
        StrucParam.Phi=30*pi/180; %Rotation angle, [rad].  
      
        %Better put the definition file in the grating profiles filefolder. 
        input_grating_file='Quadratic_sinusiodal.m';      
end

setup_grating_profile


if strcmp(StrucParam.Check_parameters,'YES')
 fprintf('The intial parameters are as below:\n\n');
 disp(StrucParam);
 
 switch StrucParam.Fourier
    case 'FFT'
        disp('Use FFT');        
   case 'Closed_form'
        disp('Use closed-form expression');        
    case 'Numerical_int'
        disp('Use Numerical integration');
end
end