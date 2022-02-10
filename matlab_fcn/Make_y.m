function [y_, y_cheb] = Make_y(Lmin,Lmax,N)
%--------------------------------------------------------------------------
% Return the physical grid points and the Chebyshev points of the domain
% Input : Range of your domain, Lmin and Lmax. 
% Output: Physical domain y_ and Chebyshev points
%--------------------------------------------------------------------------

N = N-1;
vec = (0:1:N)';
y_cheb = cos(pi*vec/N);

L_ = Lmax - Lmin;
y_ = L_*(y_cheb+1)/2;

y_ = y_ + Lmin;

end