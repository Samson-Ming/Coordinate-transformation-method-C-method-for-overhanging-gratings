%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summation of Fourier harmonics: single or multiple COSINE functions. 
% sinusodial grating is very popular;
% General representation, For some continuous grating, it profile can be approximately represented by a summation of multiple Fourier harmonics
%(Fourier series, with 0 frequency left-out, 0 frequency is merely a shift in y direction, no impact)

StrucParam.Profile='Harmonics'; % for closed-form expression only.
%
xy=textread(input_sampling_file);
xn=xy(1,:);
yn=xy(2,:);
StrucParam.dx=xn(end);
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate. 

% Rotate and resampling
Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];
xy_rot=Rot*[xn;yn];
xn=xy_rot(1,:)'; %x coordinate in the ROTATED coordinate.
yn=xy_rot(2,:)'-tan(StrucParam.Phi)*xn; %y coordinate of a(x) in the ROTATED coordinate.
N_resample=2^StrucParam.resample_order;
xq=(linspace(0,StrucParam.dx2,N_resample+1))';
yq = interp1(xn,yn,xq,'spline');

%
N_hm=length(StrucParam.Nh);

if strcmp(StrucParam.fit_method,'FFT')
    yq(1)=(yq(1)+yq(N_resample+1))/2;
    yq(N_resample+1)=[];
    F=fft(yq,N_resample)/N_resample;    % do FFT
    StrucParam.A0=F(1);
    Fn=F(2:N_hm+1); 
    StrucParam.An=2*abs(Fn);
    StrucParam.phi_n=angle(Fn);
    
elseif strcmp(StrucParam.fit_method,'Inverse')
    order=0:N_hm;
    A=cos(2*pi/StrucParam.dx2*xq*order);
    B=sin(2*pi/StrucParam.dx2*xq*order);
 %   B(:,1)=A(:,1)/2;
 %   A(:,1)=A(:,1)/2;
    B(:,1)=[];
    AB=[A,B];
    F=AB\yq;
    StrucParam.A0=F(1);
    Fn=F(2:N_hm+1)-1i*F(N_hm+2:2*N_hm+1);
    StrucParam.An=abs(Fn);
    StrucParam.phi_n=angle(Fn);    
end
%

b_x='StrucParam.A0*StrucParam.Scale+';
diff_a_x=[];
for index=1:N_hm
    b_x=[b_x,'StrucParam.An(',num2str(index),')*StrucParam.Scale*cos(2*pi*x/(StrucParam.dx2*StrucParam.Scale)*',num2str(index),'+StrucParam.phi_n(',num2str(index),'))+'];
    diff_a_x=[diff_a_x,'-StrucParam.An(',num2str(index),')*StrucParam.Scale*(2*pi/(StrucParam.dx2*StrucParam.Scale))*',num2str(index),'*sin(2*pi*x/(StrucParam.dx2*StrucParam.Scale)*',num2str(index),'+StrucParam.phi_n(',num2str(index),'))'];
end

StrucParam.b_x=b_x(1:end-1);
StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];
 %}
 