function [c, c_] = MJet(N, Lmin, Lmax)
%--------------------------------------------------------------------------
% Return energy integration matrix
% Input : Range of your domain, Lmin and Lmax. 
%         For the eigenfunction v, computes the integral of v^2 * r dr
%         within a range [0, L1]. Note that one-half is not multiplied.
% Output: Energy matrix c and Choleski decomposition c_.
% Note that norm(c_*cheb_vector)^2 yields the integral value
%--------------------------------------------------------------------------

num = round(abs(N));
c = zeros(num,num);
a = (Lmax - Lmin)/2;
b = (Lmin + Lmax)/2;

for i = (0:num-1)
    for j = (0:num-1)
        if rem(i+j,2) == 0
            p = a*b*(-((1+(-1)^(i+j))*(-1+i^2+j^2))/(i^4+(-1+j^2)^2-2*i^2*(1+j^2)));
            c(i+1,j+1)=p;
        else
            p = a*a*((-1+(-1)^(i+j))*(-4+i^2+j^2))/(i^4+(-4+j^2)^2-2*i^2*(4+j^2));
            c(i+1,j+1)=p;
        end
    end
end

c_ = chol(c);

end