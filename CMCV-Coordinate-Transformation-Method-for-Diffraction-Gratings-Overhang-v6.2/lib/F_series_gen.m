%% Generates Fouries series for given periodic function

function[f_vec] =  F_series_gen(a_fun,m,T,nDim,StrucParam) 
%input parameters - function profile, order of fft, period
%length, number of terms in the series
N=2^m;  
tol=StrucParam.accuracy; % treshold for cutting small array elements
cut_small=StrucParam.cut; % cut small array elements
y_vec=linspace(0,T,N+1); % prepare nodes
f = a_fun(y_vec);   % function values in the nodes
fhat =f;   
fhat(1)=(f(1)+f(N+1))/2;            % average of the first and last element - for faster convergence,if f(1)=f(N+1),continuous;f(1)~=f(N+1), average.
fhat(N+1)=[];

F=fft(fhat,N)/N;    % do FFT
% cut small array elements
if (cut_small==1)
    ind_small_real=(abs(real(F))<tol);
    F(ind_small_real)=1i.*imag(F(ind_small_real));
    ind_small_imag=(abs(imag(F))<tol);
    F(ind_small_imag)=real(F(ind_small_imag));
    F=CutSmallArrayElements(F,StrucParam.accuracy);
%    f_vec=f_vec(1:2:end);
end
f_vec=F(1:nDim)';    %conjagte transpose?  Becareful,here it is the conjugate, rather than the ORIGINAL Fourier coefficients.
                              %Note:f_vec(1:m2+1):conjugate of harmonics of
                              %0:m2, f_vec(m2+2:nDim):conjugate of
                              %harmonics of -m1:-1 orders.
