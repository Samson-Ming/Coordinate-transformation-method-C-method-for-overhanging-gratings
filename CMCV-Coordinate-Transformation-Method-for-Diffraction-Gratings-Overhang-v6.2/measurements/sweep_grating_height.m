H=linspace(0.3,0.45,21);
D_study_order = zeros(1,numel(H));
      if diffraction_efficiencies_c==1
        figure('Name','Reflection change with grating height'); 
        ylabel(sprintf('Efficiency of m=%2.0f reflection order', studying_order));
      elseif diffraction_efficiencies_c==2
        figure('Name','Transmission change with grating height');   
        ylabel(sprintf('Efficiency of m=%2.0f transmission order', studying_order));
      end
      xlabel('Grating height');
      box on
      hold on
      for it_index=1:numel(H)
        StrucParam.h_ORI=H(it_index);  
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
         plot(H(1:it_index),D_study_order(1:it_index));
        pause(eps);
      end
      hold off
      TotRuntime=toc;
      fprintf('Total runtime is %d sec \n',TotRuntime);