function [yPhy, yCheb] = genChebGrid(Lmin,Lmax,N)
%--------------------------------------------------------------------------
% Return the physical grid points and the Chebyshev points of the domain
% Input : Range of your domain, Lmin and Lmax. 
% Output: Physical domain y_ and Chebyshev points
%--------------------------------------------------------------------------
N = N-1;
vec = (0:1:N)';
yCheb = cos(pi*vec/N);
L = Lmax - Lmin;
yPhy = L*(yCheb+1)/2+Lmin;
end