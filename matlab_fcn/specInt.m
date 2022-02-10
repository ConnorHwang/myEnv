function integral = specInt(coeff,vecA,vecB,flag)
% -------------------------------------------------------------------------
% @ <specInt> computes volume integration of the two vectors in spectral
% space. \int r*dr dt d\theta. Note that dx is not included in this
% function.
% flag options:
% 1, the same phase with the same vectors;
% 2, different phase;
% 3, same phase but different vectors
% -------------------------------------------------------------------------
% Make sure all the inputs are in the vector format
coeff = coeff(:);
vecA = vecA(:);
vecB = vecB(:);
% Size of the vector to get N
N = size(vecA,1);
% Clenshaw-Curtis quadrature
wts = quadwts(N); wts = wts(:);
% integration
switch flag
    case 1
        temp = coeff.*abs(vecA).^2.*wts;
        temp(isinf(temp)) = 0;
        integral = 4^2*pi^2*sum(temp,'omitnan');
    case 2
        integral = 0.0;
    case 3
        temp = coeff.*(real(vecA).*real(vecB)+imag(vecA).*imag(vecB)).*wts;
        temp(isinf(temp)) = 0;
        integral = 4^2*pi^2*sum(temp,'omitnan');
end % End of switch
end