function ifft3d = ...
         spatial_ifft3DField(N1, N2, omega, n, D0_1, D0_2, evl, evc, time, axial, azim, factor)
%--------------------------------------------------------------------------
% Brief: return the 4D array, the inverse Fourier transform of
% eigenfunctions that maps in 3D plain.

% Note that when using options 3) and 4), phase adjustment has to be done.

% Input : 
% factor 1) Velocity-Vorticity Code, 2) Primitive variable code, 3) Random
% vector in fluid 1 domain, 4) Random vector in fluid 2 domain, 5)
% Interface.
% 'evc' of the form [Problem Size, nBuffer]
% 'evl' of the form [1, nBuffer]
% Output: 
%--------------------------------------------------------------------------
zi = sqrt(-1);

Nt = max(size(time,2), size(time,1));
Nz = max(size(axial,2), size(axial,1));
Nthe = max(size(azim,2), size(azim,1));
            
switch(factor)
    case 1
%         ChebToPhy = blkdiag(D0_2, D0_2, 1, D0_1, D0_1);
%         evc(2*N2+1,:) = evc(2*N2+1,:) * zi;
    case 2
        ChebToPhy = blkdiag(D0_2, D0_2, D0_2, D0_2, D0_2, D0_2, D0_2,...
                      1, 1, D0_1, D0_1, D0_1, D0_1, D0_1, D0_1, D0_1);
        evc([1:N2 4*N2+1:5*N2 7*N2+2+1:7*N2+2+N1 7*N2+2+4*N1+1:7*N2+2+5*N1 ],:)...
           = evc([1:N2 4*N2+1:5*N2 7*N2+2+1:7*N2+2+N1 7*N2+2+4*N1+1:7*N2+2+5*N1 ],:) * zi;
    case 3
        ChebToPhy = D0_1;
    case 4
        ChebToPhy = D0_2;
    otherwise
        ChebToPhy = 1;
end
    
% Pre-allocate 'ifft3d'
ifft3d = zeros(Nt, Nz, Nthe, size(evc,1));

if( n == 0 )
%--------------------------------------------------------------------------
% [ <+alp, n> <-alp, n> ]
% Should uncomment this line to use the above scheme.
    vec1 = evc(:,1);
 
        for ti = 1:Nt %parfor
            for si = 1:Nz
                for thei = 1:Nthe
                dummy1 = vec1 .* exp(2*pi*zi*evl(1)*axial(si) + 2*pi*zi*  n *azim(thei) - 2*pi*zi*omega*time(ti));    
%                 dummy1 = vec1 .* exp(zi*evl(1)*axial(si) + zi*  n *azim(thei) - zi*omega*time(ti));
                temp = dummy1 + conj(dummy1);
                ifft3d(ti,si,thei,:) = ChebToPhy * temp;
                end
            end
        end
%--------------------------------------------------------------------------
else
% -------------------------------------------------------------------------
% [ <+alp, +n> <+alp, -n> <-alp, +n> <-alp, -n> ]
test1 = evc(:,1);
test2 = evc(:,2);

for ti = 1:Nt %parfor
    for si = 1:Nz
        for thei = 1:Nthe
           dummy1 = test1 .* exp(2*pi*zi*evl(1)*axial(si) + 2*pi*zi*  n *azim(thei) - 2*pi*zi*omega*time(ti));
           dummy2 = test2 .* exp(2*pi*zi*evl(2)*axial(si) + 2*pi*zi*(-n)*azim(thei) - 2*pi*zi*omega*time(ti));
%            dummy1 = test1 .* exp(zi*evl(1)*axial(si) + zi*  n *azim(thei) - zi*omega*time(ti));
%            dummy2 = test2 .* exp(zi*evl(2)*axial(si) + zi*(-n)*azim(thei) - zi*omega*time(ti));
           temp = dummy1 + dummy2 + conj(dummy1) + conj(dummy2);
        ifft3d(ti,si,thei,:) = ChebToPhy * temp;
        end
    end
end
%--------------------------------------------------------------------------
end

end





