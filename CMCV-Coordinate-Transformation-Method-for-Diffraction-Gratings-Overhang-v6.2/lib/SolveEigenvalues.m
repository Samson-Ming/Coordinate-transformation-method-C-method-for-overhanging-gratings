function [eig1,eig2,vect1,vect2]=SolveEigenvalues(a_diff,alpha_m,beta1_m,beta2_m,StrucParam,k0)

%Normalization
alpha_m=alpha_m/k0;
beta1_m=beta1_m/k0;
beta2_m=beta2_m/k0;

a_diff_col=a_diff(StrucParam.N_Tr:1:2*StrucParam.N_Tr-1);  %positive
a_diff_row=a_diff(StrucParam.N_Tr:-1:1);       %negtive.
a_matr=toeplitz(a_diff_col,a_diff_row); % toeplitz matrix a_diff - equation (11) in the paper

Matr1=[-diag(1./(beta1_m.^2))*(diag(alpha_m)*a_matr+a_matr*diag(alpha_m)) diag(1./(beta1_m.^2))*(eye(StrucParam.N_Tr)+a_matr*a_matr);eye(StrucParam.N_Tr) zeros(StrucParam.N_Tr)];

if ~(StrucParam.n2==1i*inf)
Matr2=[-diag(1./(beta2_m.^2))*(diag(alpha_m)*a_matr+a_matr*diag(alpha_m)) diag(1./(beta2_m.^2))*(eye(StrucParam.N_Tr)+a_matr*a_matr);eye(StrucParam.N_Tr) zeros(StrucParam.N_Tr)];
end

if strcmp(StrucParam.PPT,'NO')
    
elseif strcmp(StrucParam.PPT,'YES')
    if StrucParam.decimal_digits==7
       Matr1=double(single(Matr1));  %To enhance the numerical stability with the perturbatively preconditioned C method
                              %by the k-digits preconditioning technique,
                              %original��down conversion��up conversion��original'
      if ~(StrucParam.n2==1i*inf)
      Matr2=double(single(Matr2));  %To enhance the numerical stability with the perturbatively preconditioned C method
                              %by the k-digits preconditioning technique,
                              %original��down conversion��up conversion��original
      end
    else
        
        switch StrucParam.PPT_method
            case 'LP'
        %%{
        %chop
        kk=StrucParam.decimal_digits;
        t=floor(-log2(10^-kk))+1;
        emax = 1023;
        options.round = 2;
        options.subnormal=1;
        options.format = 'c';
        options.params = [t emax];
        Matr1=chop(real(Matr1),options)+1i*chop(imag(Matr1),options);
        %}
            case 'MP'
        %%{
        mp.Digits(StrucParam.decimal_digits);
        Matr1=double(mp(Matr1));
        %}
            case 'Sym'
        %%{
        %Symbolic Math Toolbox
        digits(StrucParam.decimal_digits);
        Matr1=double(vpa(Matr1));
        %}
            case 'Round'      
        %simple round
            Matr1=round(Matr1,StrucParam.decimal_digits);  %To enhance the numerical stability with the perturbatively preconditioned C method
            otherwise
                disp('Wrong PPT method');
        end
        
       if ~(StrucParam.n2==1i*inf)
       
           switch StrucParam.PPT_method
               
            case 'LP'
            Matr2=chop(real(Matr2),options)+1i*chop(imag(Matr2),options);
            case 'MP'
            Matr2=double(mp(Matr2));
            case 'Sym'
            Matr2=double(vpa(Matr2));
            case 'Round'  
            Matr2=round(Matr2,StrucParam.decimal_digits);  %To enhance the numerical stability with the perturbatively preconditioned C method
            otherwise
                disp('Wrong PPT method');
           end
        
       end
        
    end
end


[vect1,eig1]=eig(Matr1);
eig1=1./diag(eig1);

if ~(StrucParam.n2==1i*inf)
[vect2,eig2]=eig(Matr2);
eig2=1./diag(eig2);
else
vect2=vect1*0;
eig2=([1i*ones(StrucParam.N_Tr,1);-1i*ones(StrucParam.N_Tr,1)]);
end

end