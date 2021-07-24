function [n_model]=Partial_fraction(omega,eps_inf,Omega,A)

epsilon=eps_inf;
N=numel(Omega);
for iii=1:N    
  epsilon=epsilon+A(iii)./(omega-Omega(iii))-conj(A(iii))./(omega+conj(Omega(iii)));
end

n_model=sqrt(epsilon);
end
