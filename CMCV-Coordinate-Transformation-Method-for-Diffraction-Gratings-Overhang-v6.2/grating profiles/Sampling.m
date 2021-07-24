
%%  Pure sampling points in the ORIGINAL coordinate, not fit by any function, use FFT.
StrucParam.Profile='Sampling'; % for closed-form expression only.
%StrucParam.Fourier='FFT';  %only use FFT;
%
xy=textread(input_sampling_file);
xn=xy(1,:);
yn=xy(2,:);
StrucParam.dx=xn(end);
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate. 

% Rotate and resampling
Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];
xy_rot=Rot*[xn;yn];
StrucParam.xn=xy_rot(1,:); %x coordinate in the ROTATED coordinate.
StrucParam.yn=xy_rot(2,:)-tan(StrucParam.Phi)*StrucParam.xn; %y coordinate of a(x) in the ROTATED coordinate.
%N_resample=2^StrucParam.resample_order;
%StrucParam.xn=(linspace(0,StrucParam.dx2,N_resample+1))';
%StrucParam.yn = interp1(xn,yn,StrucParam.xn,'spline');
StrucParam.dadx=(StrucParam.yn(2:end)-StrucParam.yn(1:end-1))./(StrucParam.xn(2:end)-StrucParam.xn(1:end-1))+tan(StrucParam.Phi);
StrucParam.dadx=[StrucParam.dadx,StrucParam.dadx(1)];

%%
%StrucParam.b_x='0.*(x<=StrucParam.dx2*StrucParam.Scale &x>=0)';    %Not used.
%StrucParam.diff_a_x='0.*(x<=StrucParam.dx2*StrucParam.Scale &x>=0)';    %Not used.

StrucParam.b_x='interp1(StrucParam.xn*StrucParam.Scale,StrucParam.yn*StrucParam.Scale,x,''spline'').*(x<=StrucParam.dx2*StrucParam.Scale &x>=0)';
StrucParam.diff_a_x='interp1(StrucParam.xn*StrucParam.Scale,StrucParam.dadx,x,''spline'').*(x<=StrucParam.dx2*StrucParam.Scale &x>=0)';