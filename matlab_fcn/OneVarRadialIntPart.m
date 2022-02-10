function result = OneVarRadialIntPart(N, Lmin, Lmax, Vec, rad, RadialInt, cutoff, cutofffactor)
%--------------------------------------------------------------------------
% Brief: Energy Matrix, M, to integrates the multiplication of two
% variables. Cholesky decompoisition of M will yield F.
% (M * c1_) will give the integration of c1 in physical
% space. Note, c1_ are the Chebyshev coefficient of c1.
% cutofffactor: 0) Total, 1) Upper, 2) Lower.

% Input : % RadInt - 1) Integrate in raidal direction, 2)without
%             multiplying r to the integrand.
% Output: Integrated value of a variable.
%--------------------------------------------------------------------------

[r_phy1, ~] = Make_y(Lmin,Lmax,N);
% xv = Lmin:dxv:Lmax;
xv = rad;
Ny = size(xv,2);

% Find the cutoff index
% cutoff should be within the range of 0 < cutoff < 1.
[~, minind] = min(abs(xv - cutoff));
if (1 == minind || Ny == minind)
    fprintf("in OneVarRadialIntPart function \n");
    fprintf("The cutoff location out of bound!\ncutoff should be within the range of 0 < cutoff < 1.\n");
    return
end
% Split the domain into two parts
xvDown = xv(1:minind);
xvUp = xv(minind:end);

% Interpolate integrand
FcnInterp = interp1(r_phy1, Vec, xv, 'spline');
RadInterp = interp1(r_phy1, r_phy1, xv, 'spline');

if (RadialInt == 0)
    total = trapz(xv, FcnInterp .* RadInterp);
    integration1 = trapz(xvUp, FcnInterp(minind:end) .* RadInterp(minind:end)); 
    integration2 = trapz(xvDown, FcnInterp(1:minind) .* RadInterp(1:minind)); 

else
    total = trapz(xv, FcnInterp);
    integration1 = trapz(xvUp, FcnInterp(minind:end)); 
    integration2 = trapz(xvDown, FcnInterp(1:minind)); 
end

switch(cutofffactor)
    case 0
        result = total;
    case 1
        % Upper part of the domain.
        result = integration1;
    case 2
        % Lower part of the domain.
        result = integration2;
end

end