

function [n_model]=sellmeier(lambda,A,B,C)

n_model_sq=A;
N=numel(B);
for iii=1:N
    
  n_model_sq=n_model_sq+B(iii).*lambda.^2./(lambda.^2-C(iii)^2);
end

n_model=sqrt(n_model_sq);