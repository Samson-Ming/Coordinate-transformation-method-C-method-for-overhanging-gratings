%% Divide real and imaginary eigenvalues, sort them and the respective eigenvectors
function [imag_eig1p,imag_eig2n,imag_Vec1p,imag_Vec2n]=SortEigenvaluesChand(rho1,V1,rho2,V2,tol,nDim)

%% Real
%real_eig1p_ind=(abs(imag(rho1))<tol)&(real(rho1)>tol);    %indices of eigenvalues with positive real part, incident medium
%real_eig2n_ind=(abs(imag(rho2))<tol)&(real(rho2)<-tol);    %indices of eigenvalues with negative real part, transmission medium
%real_eig1p=rho1(real_eig1p_ind);                        %positive real eigenvalues, incident medium
%real_eig2n=rho2(real_eig2n_ind);                        %negative real eigenvalues, transmission medium
%[sort_real1,idx1]=sort(real_eig1p,'descend');           %sort eigenvalues descendend, incident medium
%[sort_real2,idx2]=sort(real_eig2n,'descend');           %sort eigenvalues descendend, transmission medium
%real_eig1p=sort_real1;
%real_eig2n=sort_real2;

%% Complex

imag_eig1p_ind=(imag(rho1)>tol);                        %indices of eigenvalues with positive imaginary part, incident medium
imag_eig2n_ind=(imag(rho2)<-tol);                       %infices of eigenvalues with negative imaginary pert, transmission medium
imag_eig1p=rho1(imag_eig1p_ind);                        %eigenvalues with positive imaginary part, incident medium
imag_eig2n=rho2(imag_eig2n_ind);                        %eigenvalues with negative imaginary pert, transmission medium
[~,idx3]=sort(abs(imag(imag_eig1p)),'ascend'); %sort eigenvalues with with positive imaginary part, incident medium
[~,idx4]=sort(abs(imag(imag_eig2n)),'ascend'); %sort eigenvalues with negative imaginary pert, transmission medium
imag_eig1p=imag_eig1p(idx3);
imag_eig2n=imag_eig2n(idx4);
s_imag_Vec1p=V1(1:nDim,imag_eig1p_ind);
s_imag_Vec2n=V2(1:nDim,imag_eig2n_ind);
imag_Vec1p=s_imag_Vec1p(:,idx3);
imag_Vec2n=s_imag_Vec2n(:,idx4);