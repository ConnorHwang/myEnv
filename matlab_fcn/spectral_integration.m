function result = spectral_integration(vec1,vec2,N1,L1,coeff)
% Radially integrate the product of given vectors in spectral space.
% If the vecotr size is '1', this function assumes the integration as a 
% surface integral. This integration integrates from r = 0 to r = 1.
%--------------------------------------------------------------------------
% @ inputs
% @ vec1: input vector1
% @ vec2: input vector2
% @ L1, L2: domain extent for phase 1 and phase 2
% @ N1, N2: number of spectral points in phase 1 and phase 2
% @ coeff: coefficient of the integral
% @ Output: a complex number
%--------------------------------------------------------------------------
% make sure the vectors are columns
vec1 = vec1(:);
vec2 = vec2(:);
coeff = coeff(:);

% Check NaN/Inf value due to 1/r
vec1(isinf(vec1)) = 0;
vec2(isinf(vec2)) = 0;

% Check the vector size
sizev1 = size(vec1,1);
sizev2 = size(vec2,1);
if( sizev1 ~= sizev2 )
    fprintf('vector size must be equal!\n');
    return;
end

if( sizev1 == 1 ) % surface integral
%     result = conj(vec1).' * (coeff .* vec2);
    result = conj(real(vec1)).' * (coeff .* real(vec2)) + ...
             conj(imag(vec1)).' * (coeff .* imag(vec2));
else
    % Integration matrix
    [~,int1] = MJet(N1,0,L1);
    
    % integration
%     result = conj(int1*vec1).' * (coeff .* (int1*vec2));
    result = conj(int1*real(vec1)).' * (coeff .* (int1*real(vec2))) + ...
             conj(int1*imag(vec1)).' * (coeff .* (int1*imag(vec2)));
end
end % End of the function