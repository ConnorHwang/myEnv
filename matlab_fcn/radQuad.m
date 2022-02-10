function wi = radQuad(rgrid,N)
% Quadratures for radial integration
% rgrid: Chebyshev collocated point in radial direction
% N    : Number of points
% wi   : Weight matrix
a = rgrid(end);
b = rgrid(1);
wi = pi*(b-a)/2/(N-1)*rgrid.*(1-4/(b-a)^2*(rgrid-(a+b)/2).^2).^(1/2);
wi([1,N]) = 1e-6;
end 

% EXAMPLE TEST SCRIPT
%{
N = 100;
j = 1:N; j = j(:);
x = cos((j-1)*pi/(N-1));
a = 0;
b = 5;
r = (x+1)/2*(b-a);
f = 3*cos(r)+2*sin(r).^2+5;
figure;
plot(f,r,'*-');

phi = pi*(b-a)/2/(N-1)*r.*(1-4/(b-a)^2*(r-(a+b)/2).^2).^(1/2);
vec = sum(phi.*f);
%}