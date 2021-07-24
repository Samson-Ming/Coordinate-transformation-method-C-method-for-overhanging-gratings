function [n_model]=drude(omega,eps_inf,omega_p,gamma)

epsilon=eps_inf-omega_p^2./(omega.^2+1i*gamma*omega);

n_model=sqrt(epsilon);

end