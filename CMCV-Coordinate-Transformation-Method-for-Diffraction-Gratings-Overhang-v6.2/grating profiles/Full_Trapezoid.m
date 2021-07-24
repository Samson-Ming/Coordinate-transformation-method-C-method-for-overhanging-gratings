%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Full trapezoid, with closed_form expression.
% can also be defined as Polyline. For partial trapezoid, just use polyline profile.

StrucParam.Profile='Full_Trapezoid'; % for closed-form expression only.

StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate.

Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];

ch_rot=Rot*[StrucParam.c1_ORI,StrucParam.c2_ORI;StrucParam.h_ORI,StrucParam.h_ORI];

StrucParam.c1=ch_rot(1,1); %x coordinate in the ROTATED coordinate.
StrucParam.c2=ch_rot(1,2); %x coordinate in the ROTATED coordinate.
StrucParam.h=ch_rot(2,1)-tan(StrucParam.Phi)*StrucParam.c1; %y coordinate of a(x) in the ROTATED coordinate.

StrucParam.slope_l=StrucParam.h/StrucParam.c1;    %slope of left segment.
StrucParam.slope_r=-StrucParam.h/(StrucParam.dx2-StrucParam.c2);    %slope of left segment.

% expression of a(x') in the  the ROTATED coordinate.
% and % expression of derivitive in the  the ROTATED coordinate.
% diff_a_x=diff_b(x)/dt+tan(Phi);
% at discontinuous of derivitive, diff_a_x(x)=(diff_a_x(x-)+diff_a_x(x+))/2;
StrucParam.b_x = '(StrucParam.slope_l*x).*(x<StrucParam.c1*StrucParam.Scale & x>=0)+StrucParam.h*StrucParam.Scale.*(x<StrucParam.c2*StrucParam.Scale & x>=StrucParam.c1*StrucParam.Scale)+(StrucParam.slope_r*(x-StrucParam.dx2*StrucParam.Scale)).*(x<=StrucParam.dx2*StrucParam.Scale &x>=StrucParam.c2*StrucParam.Scale)';
StrucParam.diff_a_x = 'StrucParam.slope_l.*(x<StrucParam.c1*StrucParam.Scale & x>=0)+0.*(x<StrucParam.c2*StrucParam.Scale & x>StrucParam.c1*StrucParam.Scale)+StrucParam.slope_r.*(x<=StrucParam.dx2*StrucParam.Scale &x>StrucParam.c2*StrucParam.Scale)+(StrucParam.slope_l)/2.*(x==StrucParam.c1*StrucParam.Scale)+(StrucParam.slope_r)/2.*(x==StrucParam.c2*StrucParam.Scale)+tan(StrucParam.Phi)';