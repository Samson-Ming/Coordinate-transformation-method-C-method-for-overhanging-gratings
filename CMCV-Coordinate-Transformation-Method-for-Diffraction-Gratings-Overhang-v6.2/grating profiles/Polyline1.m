%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Polyline Type1: Profile described by Polyline in the ORIGINAL coordinate, with closed-form expression. 
%%{
% Including arbitary  triangle, arbitary trapezoid, and arbitary polygon, can be overhang.
% NOTE: an arbitary profile can also be approximated to a polyline, if only the number of sampling nodes are large enough.
% The convergence of the Polyline is more stable relative to the truncation number,  
% but the divergence problem for h/dx>>1 and h/lambda>>1 still exist, need another parametric C method.
% Both See Xihong Xu's dissertation.
% Thus it is recommended to verify your results convergence roughly by this polyline approximation for any profile.

StrucParam.Profile='Polyline'; % for closed-form expression only.

%% A polyline is generally described by the x and y coordinates of its nodes in the ORIGINAL coordinate.
 StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi); %Period in the ROTATED coordinate.
Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];

xy_rot=Rot*[StrucParam.xn_ORI;StrucParam.yn_ORI];
StrucParam.xn=xy_rot(1,:); %x coordinate in the ROTATED coordinate.
StrucParam.yn=xy_rot(2,:)-tan(StrucParam.Phi)*StrucParam.xn; %y coordinate of a(x) in the ROTATED coordinate.

StrucParam.slope=(StrucParam.yn(2:end)-StrucParam.yn(1:end-1))./(StrucParam.xn(2:end)-StrucParam.xn(1:end-1));

N_Polyline=length(StrucParam.xn)-1;

% expression of a(x') in the  the ROTATED coordinate.
% and % expression of derivitive in the  the ROTATED coordinate.
% diff_a_x=diff_b(x)/dt+tan(Phi);
% at discontinuous of derivitive, diff_a_x(x)=(diff_a_x(x-)+diff_a_x(x+))/2;

b_x=[];
diff_a_x=[];
for index=1:N_Polyline-1
    b_x=[b_x,'(StrucParam.slope(',num2str(index),')*(x-StrucParam.xn(',num2str(index),')*StrucParam.Scale)+StrucParam.yn(',num2str(index),')*StrucParam.Scale).*(x<StrucParam.xn(',num2str(index+1),')*StrucParam.Scale & x>=StrucParam.xn(',num2str(index),')*StrucParam.Scale)+'];
    diff_a_x=[diff_a_x,'StrucParam.slope(',num2str(index),').*(x<StrucParam.xn(',num2str(index+1),')*StrucParam.Scale & x>StrucParam.xn(',num2str(index),')*StrucParam.Scale)+','(StrucParam.slope(',num2str(index),')+StrucParam.slope(',num2str(index+1),'))/2.*(x==StrucParam.xn(',num2str(index+1),')*StrucParam.Scale)+'];
end

b_x=[b_x,'(StrucParam.slope(',num2str(N_Polyline),')*(x-StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)+StrucParam.yn(',num2str(N_Polyline),')*StrucParam.Scale).*(x<=StrucParam.xn(',num2str(N_Polyline+1),')*StrucParam.Scale & x>=StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)'];
diff_a_x=[diff_a_x,'StrucParam.slope(',num2str(N_Polyline),').*(x<=StrucParam.xn(',num2str(N_Polyline+1),')*StrucParam.Scale & x>StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)+','StrucParam.slope(1).*(x==StrucParam.xn(1)*StrucParam.Scale)'];

StrucParam.b_x=b_x;
StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];
