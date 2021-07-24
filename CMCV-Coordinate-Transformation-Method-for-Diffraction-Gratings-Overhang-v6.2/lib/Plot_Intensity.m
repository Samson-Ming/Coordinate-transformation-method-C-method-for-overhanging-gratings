function [ ] = Plot_Intensity(alfa0,b0,b_fun,RVec,n1_set, SB1, A, n2_set, imag_Vec1p, imag_eig1p, SB2,  imag_Vec2n,imag_eig2n, m1, m2, StrucParam)
                 
%function Plot_Intensity - generates intensity plot of the electric (TE case) or
%magnetic (TM case) field. This plot can be saved to a file. Input
%parameters are grating and incident light parameters and computed fields
%
   
    disp('Plotting intensity distribution');
    %% Plot fields
    nPoints_x=StrucParam.point_x;                            % number of nodes in x;
    d=StrucParam.dx*StrucParam.Scale;
    d2_scale=d*cos(StrucParam.Phi);
    a_fun=@ (t) b_fun(mod(t,d2_scale))+tan(StrucParam.Phi)*t; 
    %a_fun(0)
    Nx=StrucParam.Num_of_Period*nPoints_x+1;
    t=linspace(0,StrucParam.Num_of_Period*d2_scale,Nx);            % generates t-grid, i.e., the parameter t.
    y0=cos(StrucParam.Phi)*b_fun(mod(t,d2_scale));%
    x0=sec(StrucParam.Phi)*t+tan(StrucParam.Phi)*y0; 
   
    h=max(y0)-min(y0);
    nPoints_y=StrucParam.point_y;                            % number of nodes in y;
    yy=linspace(min(y0)-StrucParam.lower*h,max(y0)+StrucParam.upper*h,nPoints_y);    % generates y-grid
    xx=linspace(0,StrucParam.Num_of_Period*d,Nx);
    [X,Y]=meshgrid(xx,yy);                  % generates [X,Y] grid
    %clear yy
    
    %%
    %%  ----------------------------In the ROTATED coodinate----------------------------------------------------
    Rot=[cos(StrucParam.Phi),-sin(StrucParam.Phi);sin(StrucParam.Phi),cos(StrucParam.Phi)];
   %inv_Rot=[cos(StrucParam.Phi),sin(StrucParam.Phi);-sin(StrucParam.Phi),cos(StrucParam.Phi)];
   %[X,Y]
   X2=Rot(1,1)*X+Rot(1,2)*Y;
   Y2=Rot(2,1)*X+Rot(2,2)*Y;
    
    %wavelength=StrucParam.lambda;
    FIn=exp(1i*alfa0*X2+1i*b0*(Y2-tan(StrucParam.Phi)*X2)).*(a_fun(X2)<=Y2);  % Values of the input field on the grid
    FRPlus=zeros(nPoints_y,Nx);          %preallocate fields for positive propagation orders in superstrate medium
    FRNeg=zeros(nPoints_y,Nx);           %preallocate fields for negative propagation orders in substrate medium
    FRPlusIm=zeros(nPoints_y,Nx);        %preallocate fields for positive evanescent orders in superstrate medium
    FRNegIm=zeros(nPoints_y,Nx);         %preallocate fields for negative evanescent orders in substrate medium
    
    % computes values of reflected propagating orders on the grid 
    for mm=1:length(n1_set)
        FRPlus=(FRPlus + RVec(mm).*exp(1i*A(mm+min(n1_set)-m1).*X2 + 1i*SB1(mm+min(n1_set)-m1).*(Y2-tan(StrucParam.Phi)*X2))).*(a_fun(X2)<=Y2);
    end
     
    % computes values of reflected evanescent orders on the grid
    for mm=m1:m2
        for kk=1:length(imag_eig1p)
            FRPlusIm=FRPlusIm + (exp(1i*A(mm-m1+1).*X2).*((RVec(kk+length(n1_set))*imag_Vec1p(mm-m1+1,kk)).*exp(1i*imag_eig1p(kk).*(Y2-a_fun(X2))))).*(a_fun(X2)<=Y2);
        end
    end
    
    disp('Keep waiting, we are working on it')  % computation can be quite long
    
    % computes values of transmitted propagation orders on the grid
    for mm=1:length(n2_set)
       FRNeg=FRNeg + (RVec(mm+length(n1_set)+length(imag_eig1p)).*exp(1i*A(mm+min(n2_set)-m1).*X2 + 1i*SB2(mm+min(n2_set)-m1).*(Y2-tan(StrucParam.Phi)*X2))).*(a_fun(X2)>Y2);
    end

    % computes values of transmitted evanescent orders on the grid
    for mm=m1:m2
        for kk=1:length(imag_eig2n)
           FRNegIm=FRNegIm + (exp(1i*A(mm-m1+1).*X2).*((RVec(kk+length(n1_set)+length(imag_eig1p)+length(n2_set))*imag_Vec2n(mm-m1+1,kk)).*exp(1i*imag_eig2n(kk).*(Y2-a_fun(X2))))).*(a_fun(X2)>Y2);
        end
    end
    
    disp('It is almost done')
    
    Z=double(abs(FIn+FRPlus+FRPlusIm+FRNeg+FRNegIm));   % Computes the whole intenisty on the grid
    %Z=double(real(FIn+FRPlus+FRPlusIm+FRNeg+FRNegIm)); % Computes the values(!) of the electric/magnetic field on the grid.
    
    %% Start plot :-)
    if strcmp(StrucParam.Pol,'TM') 
    figure('Name','Distribution of intensity of |H|')
    elseif strcmp(StrucParam.Pol,'TE') 
          figure('Name','Distribution of intensity of |E|')
    end
    
    hold on
    IntensityPlot=pcolor(X/(d),Y/(d),Z);
    shading interp;
    %imagesc(xx/d,yy/d,Z);
    %IntensityPlot=contourc(X/(d),Y/(d),Z);
    colorbar;
    colormap jet;
    
    %%{
    %fy=a_fun(xx2); %plot the grating profile
    plot(x0/(d)-1,y0/(d),'w-','LineWidth',1.5); %plot the intensity
    plot(x0/(d),y0/(d),'w-','LineWidth',1.5); %plot the intensity
    xlim([x0(1)/d,x0(end)/d]);
    xlabel('x/period');
    ylabel('y/period');
    xlim([0,StrucParam.Num_of_Period]);
    ylim([yy(1),yy(end)])
    %axis equal;
    axis on
    %}
    
    %figure,plot(xx/(d),y0/(d),'k-','LineWidth',1.5); %plot the intensity
    
    % save to file if required
    if strcmp(StrucParam.save_plot,'YES')
       file_name=strcat('Intensitydistribution at wavelength=',num2str((StrucParam.lambda/StrucParam.Scale)),'.png');
       saveas(IntensityPlot, file_name); %Save plot to png file
    end
    hold off
       
end

