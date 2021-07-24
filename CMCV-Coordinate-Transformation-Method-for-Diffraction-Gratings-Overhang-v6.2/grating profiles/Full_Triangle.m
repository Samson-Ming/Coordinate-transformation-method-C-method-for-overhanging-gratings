%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Full triangle,with closed_form expression.
% can also be defined as Polyline. For partial triangle, just use polyline profile.

StrucParam.Profile='Full_Triangle'; % for closed-form expression only.
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate.

Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];

ch_rot=Rot*[StrucParam.c_ORI;StrucParam.h_ORI];
StrucParam.c=ch_rot(1); %x coordinate in the ROTATED coordinate.
StrucParam.h=ch_rot(2)-tan(StrucParam.Phi)*StrucParam.c; %y coordinate of a(x) in the ROTATED coordinate.

StrucParam.slope_l=StrucParam.h/StrucParam.c;    %slope of left segment.
StrucParam.slope_r=-StrucParam.h/(StrucParam.dx2-StrucParam.c);    %slope of left segment.

% expression of a(x') in the  the ROTATED coordinate.
% and % expression of derivitive in the  the ROTATED coordinate.
% diff_a_x=d_b(x)/dt+tan(Phi);
% at discontinuous of derivitive, diff_a_x(x)=(diff_a_x(x-)+diff_a_x(x+))/2;

StrucParam.b_x = '(StrucParam.slope_l*x).*(x<StrucParam.c*StrucParam.Scale & x>=0)+(StrucParam.slope_r*(x-StrucParam.dx2*StrucParam.Scale)).*(x<=StrucParam.dx2*StrucParam.Scale &x>=StrucParam.c*StrucParam.Scale)';
StrucParam.diff_a_x = 'StrucParam.slope_l.*(x<StrucParam.c*StrucParam.Scale & x>=0)+StrucParam.slope_r.*(x<=StrucParam.dx2*StrucParam.Scale &x>StrucParam.c*StrucParam.Scale)+(StrucParam.slope_l+StrucParam.slope_r)/2.*(x==StrucParam.c*StrucParam.Scale)+tan(StrucParam.Phi)';
%}

%% Below are some examples.
% Put them in the SetInitialParameters.m

%% Example3.2: Full isosceles triangle
%{
StrucParam.dx = 2*wavelength;
StrucParam.Phi=0*pi/180;

StrucParam.c=StrucParam.dx/2;  % x coordinate of the vertex;
StrucParam.h=tand(30)*StrucParam.dx/2;
%}

%% Example3.3: Full right triangle
%{
StrucParam.dx = 2*StrucParam.wavelength;
StrucParam.Phi=acot(2*tand(30));

StrucParam.c=StrucParam.dx;  % x coordinate of the vertex;
StrucParam.h=tand(30)*StrucParam.dx;
%}