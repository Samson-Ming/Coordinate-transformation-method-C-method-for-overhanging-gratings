%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PolylineType2: Profile described by Polyline in the ROTATED coordinate, with closed-form expression. 
%%{
% % Including any NON-overhang arbitary  triangle, arbitary trapezoid, and arbitary polygon. NOTE still be a function.
% NOTE: an arbitary profile (a function) can also be approximated to a polyline, if only the number of sampling nodes are large enough.
% The convergence of the Polyline is more stable relative to the truncation number,  
% but the divergence problem for h/dx>>1 and h/lambda>>1 still exist, need another parametric C method.
% Both See Xihong Xu's dissertation.
% Thus it is recommended to verify your results convergence roughly by this polyline approximation for any profile, see the sinusodial grating example.

StrucParam.Profile='Polyline'; % for closed-form expression only.

%% Here the polyline is described by the x' and y' coordinates of the nodes of a(x') in the ROTATED coordinate.

StrucParam.slope=(StrucParam.yn(2:end)-StrucParam.yn(1:end-1))./(StrucParam.xn(2:end)-StrucParam.xn(1:end-1));

N_Polyline=length(StrucParam.xn)-1;
b_x=[];
diff_a_x=[];
for index=1:N_Polyline-1
    b_x=[b_x,'(StrucParam.slope(',num2str(index),')*(x-StrucParam.xn(',num2str(index),')*StrucParam.Scale)+StrucParam.yn(',num2str(index),')*StrucParam.Scale).*(x<StrucParam.xn(',num2str(index+1),')*StrucParam.Scale & x>=StrucParam.xn(',num2str(index),')*StrucParam.Scale)+'];
    diff_a_x=[diff_a_x,'StrucParam.slope(',num2str(index),').*(x<StrucParam.xn(',num2str(index+1),')*StrucParam.Scale & x>StrucParam.xn(',num2str(index),')*StrucParam.Scale)+','(StrucParam.slope(',num2str(index),')+StrucParam.slope(',num2str(index+1),'))/2.*(x==StrucParam.xn(',num2str(index+1),')*StrucParam.Scale)+'];
end

% expression of a(x') in the  the ROTATED coordinate.
% and % expression of derivitive in the  the ROTATED coordinate.
% diff_a_x=diff_b(x)/dt+tan(Phi);
% at discontinuous of derivitive, diff_a_x(x)=(diff_a_x(x-)+diff_a_x(x+))/2;

b_x=[b_x,'(StrucParam.slope(',num2str(N_Polyline),')*(x-StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)+StrucParam.yn(',num2str(N_Polyline),')*StrucParam.Scale).*(x<=StrucParam.xn(',num2str(N_Polyline+1),')*StrucParam.Scale & x>=StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)'];
diff_a_x=[diff_a_x,'StrucParam.slope(',num2str(N_Polyline),').*(x<=StrucParam.xn(',num2str(N_Polyline+1),')*StrucParam.Scale & x>StrucParam.xn(',num2str(N_Polyline),')*StrucParam.Scale)+','StrucParam.slope(1).*(x==StrucParam.xn(1)*StrucParam.Scale)'];

StrucParam.b_x=b_x;
StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];
%}


%% Below are some examples.
% Put them in the SetInitialParameters.m

% Rotated triangle
%{
StrucParam.dx = 0.425;
StrucParam.Phi=acot(0.36/(0.4211-0.425/2));
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi);

 StrucParam.xn=[0   0.053471815577109   0.183863443529703   0.314255071482297   StrucParam.dx2]; % The 1st must be 0, the last must be period, and xn(m)<xn(m+1);
 StrucParam.yn=[0 0 0.36 0 0];
%}