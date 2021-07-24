function [n_model]=lorentz(omega,eps_inf,omega_p,omega0,gamma)

epsilon=eps_inf;
N=numel(omega_p);
for iii=1:N    
  epsilon=epsilon+omega_p(iii)^2./(omega.^2-omega0(iii)^2+1i*gamma(iii)*omega);
end

n_model=sqrt(epsilon);
end
