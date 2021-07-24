%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summation of Fourier harmonics: single or multiple COSINE functions. 
% sinusodial grating is very popular;
% General representation, For some continuous grating, it profile can be approximately represented by a summation of multiple Fourier harmonics
%(Fourier series, with 0 frequency left-out, 0 frequency is merely a shift in y direction, no impact)

StrucParam.Profile='Harmonics'; % for closed-form expression only.
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate.     

N_hm=length(StrucParam.An);
b_x='StrucParam.A0*StrucParam.Scale+';
diff_a_x=[];
for index=1:N_hm
    b_x=[b_x,'StrucParam.An(',num2str(index),')*StrucParam.Scale*cos(2*pi*x/(StrucParam.dx2*StrucParam.Scale)*',num2str(index),'+StrucParam.phi_n(',num2str(index),'))+'];
    diff_a_x=[diff_a_x,'-StrucParam.An(',num2str(index),')*StrucParam.Scale*(2*pi/(StrucParam.dx2*StrucParam.Scale))*',num2str(index),'*sin(2*pi*x/(StrucParam.dx2*StrucParam.Scale)*',num2str(index),'+StrucParam.phi_n(',num2str(index),'))'];
end

StrucParam.b_x=b_x(1:end-1);
StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];
 %}
 
 %% Below are some examples.
% Put them in the SetInitialParameters.m

% Rotated triangle
%{
%StrucParam.An=0.360/2/cos(StrucParam.Phi);  % Fourier-coefficients amplitude
%StrucParam.phi_n=0*pi; % Fourier-coefficients phase%
%}