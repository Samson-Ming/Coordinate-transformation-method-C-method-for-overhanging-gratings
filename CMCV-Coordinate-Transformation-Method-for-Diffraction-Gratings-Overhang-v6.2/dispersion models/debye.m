function [n_model]=debye(omega,eps_inf,delta,tau)

epsilon=eps_inf;
N=numel(delta);
for iii=1:N
epsilon=epsilon+delta(iii)./(1-1i*omega*tau(iii));
end

n_model=sqrt(epsilon);

end