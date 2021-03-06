function [n_model]=rix_spline(lambda,material)

%filename=[material];


all_path=[pwd,'/dispersion models/',material];

data_s=dlmread(all_path);

velikost=size(data_s);

% kontrola zadání

% ? is lambda(1)<lambda(2)
if data_s(1,1)>data_s(2,1)
% -> revert
data_s=flipud(data_s);
end

if min(lambda)>=data_s(1,1) && max(lambda)<=data_s(velikost(1),1)
else
    error('lambda is out of range')
end
%-------------------------------------------------------------------------

data_lambda=data_s(:,1);
real_rix=data_s(:,2);

if size(data_s,2)==2
    imag_rix=zeros(size(data_s,1),1);
else
    imag_rix=data_s(:,3);
end

n_r=spline(data_lambda,real_rix,lambda);
n_i=spline(data_lambda,imag_rix,lambda);

n_model=n_r+n_i*1i;