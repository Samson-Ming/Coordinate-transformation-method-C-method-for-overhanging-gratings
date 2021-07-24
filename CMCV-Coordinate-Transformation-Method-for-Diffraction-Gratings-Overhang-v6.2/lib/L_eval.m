function L=L_eval(gamma, m, StrucParam)

if strcmp(StrucParam.Fourier,'Closed_form')
%% Closed-form integration
%Tip: you can use the symbollic compution toolbox in MATLAB or any other software to help obtain 
%the analytical expression of FF quikley and accurately first, and then copy it into codes here.

switch StrucParam.Profile
    case 'Polyline'
          
        K=2*pi/(StrucParam.dx2*StrucParam.Scale);
                
        L_segment=1/(StrucParam.dx2*StrucParam.Scale)*(exp(1i*gamma*(StrucParam.yn(1:end-1)*StrucParam.Scale-StrucParam.slope.*StrucParam.xn(1:end-1)*StrucParam.Scale)).*...
        1./(1i*(gamma*StrucParam.slope-m*K)).*...
        (exp(1i*(gamma*StrucParam.slope-m*K).*StrucParam.xn(2:end)*StrucParam.Scale)-exp(1i*(gamma*StrucParam.slope-m*K).*StrucParam.xn(1:end-1)*StrucParam.Scale)));      
            
         index=find(gamma*StrucParam.slope-m*K==0);
          
         if ~isempty(index)
            for iii=1:length(index)
                L_segment(index(iii))=1/(StrucParam.dx2*StrucParam.Scale).*(exp(1i*gamma*(StrucParam.yn(index(iii))*StrucParam.Scale-StrucParam.slope(index(iii)).*StrucParam.xn(index(iii))*StrucParam.Scale)))....
                                      .*(StrucParam.xn(index(iii)+1)-StrucParam.xn(index(iii)))*StrucParam.Scale;
            end
         end
         
         L=sum(L_segment);                    
        
    case 'Full_Triangle'        
        L=(exp(1i*(gamma*StrucParam.slope_l-2*m*pi/(StrucParam.dx2*StrucParam.Scale))*StrucParam.c*StrucParam.Scale)-1)/(1i*(StrucParam.dx2*StrucParam.Scale)*(gamma*StrucParam.slope_l-2*m*pi/(StrucParam.dx2*StrucParam.Scale)))...
          +exp(i*gamma*(StrucParam.dx2*StrucParam.Scale)*(-StrucParam.slope_r))*(exp(-i*gamma*(-StrucParam.slope_r)*(StrucParam.dx2*StrucParam.Scale))-exp(-i*(gamma*(-StrucParam.slope_r)+2*m*pi/(StrucParam.dx2*StrucParam.Scale))*StrucParam.c*StrucParam.Scale))...
          ./(-i*(StrucParam.dx2*StrucParam.Scale)*(gamma*(-StrucParam.slope_r)+2*m*pi/(StrucParam.dx2*StrucParam.Scale)));
      
    case 'Full_Trapezoid'
        K=2*pi/(StrucParam.dx2*StrucParam.Scale);
        if m==0
        L=1/(StrucParam.dx2*StrucParam.Scale)*(...
            1./(1i*(gamma*StrucParam.slope_l)).*(exp(1i*(gamma*StrucParam.slope_l).*StrucParam.c1*StrucParam.Scale)-1)+...
            exp(1i*gamma*StrucParam.h*StrucParam.Scale).*(StrucParam.c2-StrucParam.c1)*StrucParam.Scale+...
            exp(1i*gamma*(StrucParam.h*StrucParam.Scale-StrucParam.slope_r.*StrucParam.c2*StrucParam.Scale)).*1./(1i*(gamma*StrucParam.slope_r)).*...
        (exp(1i*(gamma*StrucParam.slope_r).*StrucParam.dx2*StrucParam.Scale)-exp(1i*(gamma*StrucParam.slope_r).*StrucParam.c2*StrucParam.Scale))...
         );
        else
        L=1/(StrucParam.dx2*StrucParam.Scale)*(...
            1./(1i*(gamma*StrucParam.slope_l-m*K)).*(exp(1i*(gamma*StrucParam.slope_l-m*K).*StrucParam.c1*StrucParam.Scale)-1)+...
            exp(1i*gamma*(StrucParam.h*StrucParam.Scale)).*1./(1i*(-m*K)).*(exp(1i*(-m*K).*StrucParam.c2*StrucParam.Scale)-exp(1i*(-m*K).*StrucParam.c1*StrucParam.Scale))+...                                      
            exp(1i*gamma*(StrucParam.h*StrucParam.Scale-StrucParam.slope_r.*StrucParam.c2*StrucParam.Scale)).*1./(1i*(gamma*StrucParam.slope_r-m*K)).*...
           (exp(1i*(gamma*StrucParam.slope_r).*StrucParam.dx2*StrucParam.Scale)-exp(1i*(gamma*StrucParam.slope_r-m*K).*StrucParam.c2*StrucParam.Scale))...
           );
        end
      
        
    case 'Harmonics'
              N_hm=length(StrucParam.An);
         if N_hm>1
            StrucParam.Fourier='Numerical_int';
         else
           L=exp(1i*gamma*StrucParam.A0*StrucParam.Scale)*1i^m*besselj(m,gamma*StrucParam.An*StrucParam.Scale)*exp(1i*m*StrucParam.phi_n); 
         end
         
    otherwise
        disp('Currently no Closed_form expression for this profile, switch to the Numerical_int');
        StrucParam.Fourier='Numerical_int';
end

end

if strcmp(StrucParam.Fourier,'Numerical_int')
%% Numerical integration

FF=@(x) (1/(StrucParam.dx2*StrucParam.Scale))*exp(1i*gamma*eval(StrucParam.b_x) - 1i*m*2*pi*x/(StrucParam.dx2*StrucParam.Scale));
L=quadgk(FF, 0, StrucParam.dx2*StrucParam.Scale);

end

end

