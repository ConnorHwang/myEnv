function [RadialIntField, Integration] = ifft3DIntPart_original(ifft3d, time, axial, azim, N, Lmin, Lmax, dxv, IntDir, RadInt, cutoff, cutofffactor)
%--------------------------------------------------------------------------
% Brief: 
%
% Input: 4D Matrix of which the form (time, axial, azimuthal, radial)
% IntDir - 0) Do nothing, 1)time, 2)axial, 3)azimuthal
% RadInt - 0) Do not integrate, 1) Integrate in raidal direction, 2)without
%             multiplying r to the integrand.
% Output: 
% RadialIntField - 3D Matrix in radial direction
% Integration - 3D Matrix integrated either in time or in space
%--------------------------------------------------------------------------
Nt   = size(ifft3d, 1);
Nz   = size(ifft3d, 2);
Nthe = size(ifft3d, 3);

RadialIntField = zeros(Nt,Nz,Nthe);

if (RadInt == 1)
    for i = 1:Nt % parfor
        for j = 1:Nz
            for k = 1:Nthe
            temp = squeeze(ifft3d(i,j,k,:));
            RadialIntField(i,j,k) = OneVarRadialIntPart(N, Lmin, Lmax, temp, dxv, 0, cutoff, cutofffactor); 
            end
        end
    end
elseif (RadInt == 2)
    for i = 1:Nt % parfor
        for j = 1:Nz
            for k = 1:Nthe
            temp = squeeze(ifft3d(i,j,k,:));
            RadialIntField(i,j,k) = OneVarRadialIntPart(N, Lmin, Lmax, temp, dxv, 1, cutoff, cutofffactor); 
            end
        end
    end
else
    % If we don't have to integrate over the radial direction
    RadialIntField = ifft3d;
end

switch(IntDir)
    case 0
        % Do nothing
        Integration = -1;
    case 1
%         dt = time(2) - time(1);
%         tBuffer = sum((RadialIntField(1:Nt-1,:,:) + RadialIntField(2:Nt,:,:))/2 * dt,1);
        tBuffer = trapz(time, RadialIntField, 1);
        Integration = tBuffer;
    case 2
%         dz = axial(2) - axial(1);
%         zBuffer = sum((RadialIntField(:,1:Nz-1,:) + RadialIntField(:,2:Nz,:))/2 * dz,2);
        zBuffer = trapz(axial, RadialIntField, 2);
        Integration = zBuffer;
    otherwise
%         dthe = azim(2)- azim(1);
%         theBuffer = sum((RadialIntField(:,:,1:Nthe-1) + RadialIntField(:,:,2:Nthe))/2 * dthe,3);
        theBuffer = trapz(azim, RadialIntField, 3);
        Integration = theBuffer;
end
end % End of the function