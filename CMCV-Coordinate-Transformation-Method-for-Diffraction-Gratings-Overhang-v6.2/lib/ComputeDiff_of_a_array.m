function [a_diff]=ComputeDiff_of_a_array(StrucParam)

if strcmp(StrucParam.Fourier,'Closed_form')
%% Closed-form
%Tip: you can use the symbollic compution toolbox in MATLAB or any other software to help obtain 
%the analytical expression of FF quikley and accurately first, and then copy it into codes here.
a_diff=zeros(1,StrucParam.N_Tr);

switch StrucParam.Profile
    case 'Polyline'
        it=1;
        %for n=n_set(1):1:n_set(end)
        for n=0:StrucParam.N_Tr-1
            if n==0
               a_diff(it)=1/StrucParam.dx2*sum(StrucParam.slope.*(StrucParam.xn(2:end)-StrucParam.xn(1:end-1))); 
            else
                a_diff(it)=1i/(2*n*pi)*sum(StrucParam.slope.*(exp(-1i*n*2*pi/StrucParam.dx2.*StrucParam.xn(2:end))-exp(-1i*n*2*pi/StrucParam.dx2.*StrucParam.xn(1:end-1))));
            end
            it=it+1;  
        end
        a_diff(1)=a_diff(1)+tan(StrucParam.Phi);
  
    case 'Full_Triangle'
        it=1;
        for n=0:StrucParam.N_Tr-1
            if n==0       
            a_diff(it)=StrucParam.c/StrucParam.dx2*(StrucParam.slope_l-StrucParam.slope_r)+StrucParam.slope_r;
            else         
            a_diff(it)=1i/(2*n*pi)*(StrucParam.slope_l-StrucParam.slope_r)*(exp(-1i*n*2*pi*StrucParam.c/StrucParam.dx2)-1);
            end
            it=it+1;  
        end
        a_diff(1)=a_diff(1)+tan(StrucParam.Phi);
        
    case 'Full_Trapezoid'       
        it=1;
        for n=0:StrucParam.N_Tr-1
            if n==0
               a_diff(it)=1/StrucParam.dx2*(StrucParam.slope_l*StrucParam.c1+StrucParam.slope_r*(StrucParam.dx2-StrucParam.c2));
            else                           
                a_diff(it)=1i/(2*n*pi)*...
                          (StrucParam.slope_l.*(exp(-1i*n*2*pi/StrucParam.dx2.*StrucParam.c1)-1)+...
                          StrucParam.slope_r.*(1-exp(-1i*n*2*pi/StrucParam.dx2.*StrucParam.c2)));
            end
            it=it+1;  
        end
        a_diff(1)=a_diff(1)+tan(StrucParam.Phi);
        
        
    case 'Harmonics'      
         N_hm=length(StrucParam.An);
         it=1;
         for n=0:StrucParam.N_Tr-1
             for idx=1:N_hm
                 if n==idx
                     a_diff(it)=-StrucParam.An(idx)*StrucParam.Scale*idx*2*pi/(StrucParam.dx2*StrucParam.Scale)*exp(1i*StrucParam.phi_n(idx))/(2i);
            %     elseif n==-idx
            %        a_diff(it)=StrucParam.An(idx)*StrucParam.Scale*idx*2*pi/(StrucParam.dx2*StrucParam.Scale)*exp(-1i*StrucParam.phi_n(idx))/(2i);
                 end
             end
             it=it+1;
         end
         a_diff(1)=a_diff(1)+tan(StrucParam.Phi);

    otherwise
        disp('Currently no Closed_form expression for this profile, switch to the Numerical_int');
        StrucParam.Fourier='Numerical_int';
end

a_diff=[conj(fliplr(a_diff)),a_diff(2:end)];

end

if strcmp(StrucParam.Fourier,'Numerical_int')
%% Numerical integration
a_diff=zeros(1,StrucParam.N_Tr);
it=1;
%for n=n_set(1):1:n_set(end)
for n=0:StrucParam.N_Tr-1
    F=@(x)(1/(StrucParam.dx2*StrucParam.Scale))*exp(-1i*n*2*pi*x/(StrucParam.dx2*StrucParam.Scale)).*eval(StrucParam.diff_a_x);
    a_diff(it)=quadgk(F,0,(StrucParam.dx2*StrucParam.Scale));
    %x=linspace(0,StrucParam.dx2,1e4);
    %y=F(x);
    %a_diff(it)=trapz(x,y);
    it=it+1;
end
   a_diff=[conj(fliplr(a_diff)),a_diff(2:end)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if StrucParam.cut==1 %cut small elements in the array
    a_diff=CutSmallArrayElements(a_diff,StrucParam.accuracy);
end

end