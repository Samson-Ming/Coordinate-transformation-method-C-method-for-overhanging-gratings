function out=Hz2nm(in)
c0 = 2.99792458e10; % speed of light  
nm = 1e-7;          % nanometers in CGS
out=c0./in./nm;
end