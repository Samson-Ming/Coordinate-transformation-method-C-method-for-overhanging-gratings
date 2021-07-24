%a(x_rot)=(h2.*[1-(1-sin(pi/d2.*x_rot).^2alpha).^beta];

StrucParam.Profile='Power_sine'; % for closed-form expression only.
StrucParam.dx2=StrucParam.dx.*cos(StrucParam.Phi); %Period in the ROTATED coordinate.
StrucParam.h2=StrucParam.h/cos(StrucParam.Phi); %Height in the ROTATED coordinate.

        
StrucParam.b_x='(StrucParam.h2.*StrucParam.Scale).*(1-(1-sin(pi/(StrucParam.dx2.*StrucParam.Scale).*x).^(2.*StrucParam.Alpha)).^StrucParam.Beta)';

diff_a_x=['-(StrucParam.h2.*StrucParam.Scale).*StrucParam.Beta.*(1-sin(pi/(StrucParam.dx2.*StrucParam.Scale).*x).^(2.*StrucParam.Alpha)).^(StrucParam.Beta-1)',...
              '.*(-2.*StrucParam.Alpha).*sin(pi/(StrucParam.dx2.*StrucParam.Scale).*x).^(2.*StrucParam.Alpha-1)',...
              '.*cos(pi/(StrucParam.dx2.*StrucParam.Scale).*x).*pi/(StrucParam.dx2.*StrucParam.Scale)'];
StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];
         