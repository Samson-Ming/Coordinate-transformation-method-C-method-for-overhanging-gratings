 %% Other user-defined profiles defined by any continuous differential functions.
 %% BE CAREFUL with the StrucParam.Scale term.
 % 0. For theoritical study purpose, you can test any continuous differential functions as you like.
 % 1.For real grating profile, if you have some fine sampling points of the profile, you can just use the Polyline (grating=0) profile or by its Fourier components (grating=5) or pure sampling dots (grating=6).
 % 3.You can also use some funtions like B-spline,cubic-spline,Legendre/Chebyshev/Gegenbauer polynomial expansion to fit the profile.

 %% Example:Quadratic sinusiodal profile
StrucParam.dx2=StrucParam.dx*cos(StrucParam.Phi);
StrucParam.h=0.2008; %depth of the grating in the ROTATED coordinate.
% Shape of the grating
h2=1.083;  % unit:radius                    
StrucParam.aa=(pi-h2)*(2/StrucParam.dx2); %unit:1/length, then multiplily 1/StrucParam.Scale
StrucParam.bb=h2*(2/StrucParam.dx2)^2;  % unit:1/length^2, then multiplily 1/StrucParam.Scale^2
StrucParam.b_x=[' (StrucParam.h*StrucParam.Scale/2)*cos(StrucParam.aa/StrucParam.Scale.*x+StrucParam.bb/StrucParam.Scale^2.*(x.^2)).*(x<=StrucParam.dx2*StrucParam.Scale/2)+(StrucParam.h*StrucParam.Scale/2)*cos(StrucParam.aa/StrucParam.Scale.*(-x+',...
                         'StrucParam.dx2*StrucParam.Scale)+StrucParam.bb/StrucParam.Scale^2.*(-x+StrucParam.dx2*StrucParam.Scale).^2).*(x>=StrucParam.dx2*StrucParam.Scale/2)'];
diff_a_x=['-(StrucParam.h*StrucParam.Scale/2)*(StrucParam.aa/StrucParam.Scale+2*StrucParam.bb/StrucParam.Scale^2.*x).*sin(StrucParam.aa/StrucParam.Scale.*x+StrucParam.bb/StrucParam.Scale^2.*(x.^2)).*(x<=StrucParam.dx2*StrucParam.Scale/2)+',...
              '(StrucParam.h*StrucParam.Scale/2)*(StrucParam.aa/StrucParam.Scale+2*StrucParam.bb/StrucParam.Scale^2.*(-x+StrucParam.dx2*StrucParam.Scale)).*sin(StrucParam.aa/StrucParam.Scale.*(-x+StrucParam.dx2*StrucParam.Scale)+',...
              'StrucParam.bb/StrucParam.Scale^2.*(-x+StrucParam.dx2*StrucParam.Scale).^2).*(x>=StrucParam.dx2*StrucParam.Scale/2)'];

StrucParam.diff_a_x=[diff_a_x,'+tan(StrucParam.Phi)'];