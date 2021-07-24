switch StrucParam.measurement
    case 0
       %% Single point Diffraction efficiency        
        [R_tot,T_tot,R_orders,T_orders,R,T,m_set,n1_set_ind,n2_set_ind,beta1_m_rot,beta2_m_rot,beta1_m_anti_rot,beta2_m_anti_rot,k0,eig1_p,  eig1_m,eig2_p,  eig2_m] = C_method_compute(StrucParam);    %R_orders,T_orders,
        timerVal = toc;
        fprintf('Reflection efficiency for %s polarization\n', StrucParam.Pol);
        if strcmp(StrucParam.Pol,'BOTH')
            fprintf(' First colume: orders; Second colume:reflection efficiency for TM; Third colume:reflection efficiency for TE\n');
        else
            fprintf(' First colume: orders; Second colume:reflection efficiency\n');
        end
        disp(R_orders)
        if strcmp(StrucParam.Pol,'BOTH')
            fprintf(' First colume: total reflection efficiency for TM; Second colume: total reflection efficiency for TE\n');
        end
        disp('Total reflection efficiency');
        disp(R_tot)
        fprintf('Transmission efficiency for %s polarization\n', StrucParam.Pol);
         if strcmp(StrucParam.Pol,'BOTH')
            fprintf(' First colume: orders; Second colume:transmission efficiency for TM; Third colume:transmission efficiency for TE\n');
         else
            fprintf(' First colume: orders; Second colume:transmission efficiency\n');
        end
        disp(T_orders)
        if strcmp(StrucParam.Pol,'BOTH')
            fprintf(' First colume: total transmission efficiency for TM; Second colume: total transmission efficiency for TE\n');
        end
        disp('Total transmission efficiency');
        disp(T_tot)
        
        if strcmp(StrucParam.Pol,'BOTH')
            fprintf(' First colume: total efficiency for TM; Second colume: total efficiency for TE\n');
        end
        disp('R_tot+T_tot');
        disp(R_tot+T_tot)
        
        fprintf('Runtime is %d sec \n\n',timerVal);
        fprintf('The real truncted number is %d . \n\n', StrucParam.N_Tr);
        
    case 1
    %% Convergence test, sweep N_Tr
    Convergence = zeros(1,numel(N));
   figure('Name','Convergence test');
       xlabel('Real trunction number');
    if diffraction_efficiencies_c==1
        ylabel(sprintf('Convergence of efficiency of m=%2.0f reflection order', studying_order));
    elseif diffraction_efficiencies_c==2
        ylabel(sprintf('Convergence of efficiency of m=%2.0f transmission order', studying_order));
    end
   hold on
   box on
   R_convergence=[];
   T_convergence=[];
   h_waitbar=waitbar(0,'Convergence test');
for it_index=1:numel(N)
    StrucParam.N_Tr = 2*N(it_index) + 1;
     [R_tot,T_tot,R_orders,T_orders,R,T,m_set,n1_set_ind,n2_set_ind] = C_method_compute(StrucParam);
     if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
     end
    if diffraction_efficiencies_c==1
        Convergence(it_index)=R(studying_order+1-m_set(1),Pol_ind);        
    elseif diffraction_efficiencies_c==2
        Convergence(it_index)=T(studying_order+1-m_set(1),Pol_ind);        
    end
    R_convergence=[R_convergence,R(n1_set_ind,Pol_ind)];
    if ~isempty(n2_set_ind)
        T_convergence=[T_convergence,T(n2_set_ind,Pol_ind)];
    end
    %T_convergence=[T_convergence,T(n2_set_ind,Pol_ind)];
    plot(2*N(1:it_index) + 1,(Convergence(1:it_index)));
    waitbar(it_index/numel(N),h_waitbar);
    pause(eps);
end
hold off
TotRuntime=toc;
fprintf('Total runtime is %d sec \n',TotRuntime);
   
    case 2
      %% Sweep wavelength.
      %% Update StrucParam.Scale if wavelength or wavenumber is normalized.
      %% Update StrucParam if dispersive.
        %Dispersion: n1=n1(Wavelength); n2=n2(Wavelength); StrucParam = SetConstantsByPolarization(StrucParam,n1, n2, Polarization);
        %NOTE that the wavelength (or frequency) is non-scaled in despersion model.
        %Material dispersion is defined in the "dispersion model" filefolder, but currently not avaible, will refine it in later version.
        D_study_order = zeros(1,numel(Wavelength));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection spectra'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission spectra');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('\lambda');
      hold on
      box on
      N=numel(n2);
      for it_index=1:numel(Wavelength)
        %Dispersion: n1=n1(Wavelength); n2=n2(Wavelength); StrucParam = SetConstantsByPolarization(StrucParam,n1, n2, Polarization);  
        StrucParam.Scale=1/Wavelength(it_index); %Nomalize the wavelength to 1 (no unit).
        StrucParam.lambda =Wavelength(it_index)*StrucParam.Scale; %wavelength need scaling here.
        
        %dispersion
        if N>1
            StrucParam = SetConstantsByPolarization(StrucParam,n1, n2(it_index),Polarization); 
        end
            
        [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
        if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
        end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(Wavelength(1:it_index),D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime);   
            
    case 3
        % 3 --- dependence of the diffraction efficiencies on the incident angle theta
        D_study_order = zeros(1,numel(Theta));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with incident angle'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with incident angle');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('\theta (degree)');
      hold on
      box on
      
       for it_index=1:numel(Theta)
            StrucParam.theta=Theta(it_index);
           [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
           if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
           end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(Theta(1:it_index)*180/pi,D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime); 
                  
    case 4
        %% Update StrucParam.
         %dependence of the diffraction efficiencies on the refractive index $n_1$ ($n_g$). Update StrucParam if dispersive.
         D_study_order = zeros(1,numel(n1));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with refraction index of the first media'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with refraction index of the first media');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('n_1');
      hold on
      box on
      
       for it_index=1:numel(n1)
        StrucParam = SetConstantsByPolarization(StrucParam,n1(it_index), n2, Polarization);  
        [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
        if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
        end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(n1(1:it_index),D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime);   
        
    case 5
        %% Update StrucParam.
         %dependence of the diffraction efficiencies on the refractive index $n_2$ ($n_r$). Update StrucParam if dispersive.
         D_study_order = zeros(1,numel(n2));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with refraction index of the second media'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with refraction index of the second media');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('n_2');
      hold on
      box on
      
       for it_index=1:numel(n2)
        StrucParam = SetConstantsByPolarization(StrucParam,n1, n2(it_index), Polarization);  
        [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
        if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
        end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(n2(1:it_index),D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime);   
        
    case 6
          %% Inapplicable for grating=0,5,6,and maybe 7.
        %% Update StrucParam.Scale if period is normalized.
        D_study_order = zeros(1,numel(Period));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with period'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with period');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('Period');
      hold on
      box on
      
      for it_index=1:numel( Period)
        %StrucParam.Scale=1/Period(it_index); %Nomalize the wavelength to 1 (no unit).
        %StrucParam.lambda =StrucParam.wavelength*StrucParam.Scale; %wavelength need scaling here.
        StrucParam.dx=Period(it_index);  %DO NOT scale length variable except for wavelength here.
        setup_grating_profile;  %Need to update the grating profile.
        [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
        if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
        end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(Period(1:it_index),D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime);   index
        
    case 7
        %% ONLY some profile disribed in the ORIGINAL coordinate (grating=0,2,3,5,6, and maybe 7)
        %% There might be problem if the a(x) is not function after rotation.
        %% If Phi is in the correct range, this curve should be a flat line 
        D_study_order = zeros(1,numel(Phi));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with rotation angle'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with rotation angle');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('Rotation angle \Phi (degree)');
      hold on
      box on
      
       for it_index=1:numel(Phi)
           %SetInitialParameters
            StrucParam.Phi=Phi(it_index);
            setup_grating_profile;  %Need to update the grating profile.
            StrucParam.Scale=1/StrucParam.wavelength; %Nomalize the wavelength to 1 (no unit).
            StrucParam.lambda =StrucParam.wavelength*StrucParam.Scale; %wavelength need scaling here.
           [R_tot,T_tot,R_orders,T_orders,R,T,m_set] = C_method_compute(StrucParam);
           if ~ismember(studying_order,m_set)
         studying_order=0;
         disp('study order out of truncated range,set to 0th order');
         title('Warning: study order has already been changed to the 0th order');
           end
        
         if diffraction_efficiencies_c==1
            D_study_order (it_index)=R(studying_order+1-m_set(1),Pol_ind);
            elseif diffraction_efficiencies_c==2
            D_study_order (it_index)=T(studying_order+1-m_set(1),Pol_ind);
         end  
         plot(Phi(1:it_index)*180/pi,D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime); 
        
    case 8
        %% DO NOT use "index" as iteration indicator because it has alread been used in SetInitialParameters.m
        eval(measurement_file(1:end-2));
        
end