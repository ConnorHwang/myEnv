function [A,B,C,E] = genOpGlobal(Lx,Nx,Ly,Ny,xgrid,delta)
global kz kx zi Re U dyU dyyU V dyV
% generateOperators
[~,D0x,D1x,D2x,~,~] = genChebGlobal(Nx,Lx);
% [~,D0x] = fourdif(Nx,0);
% [~,D1x] = fourdif(Nx,1);
% [~,D2x] = fourdif(Nx,2);
[~,D0y,D1y,D2y,~,~] = genChebGlobal(Ny,Ly);
%--------------------------------------------------------------------------
% Parallel test
if(Nx == 1)
    D1x = zi*kx;
    D2x = -kx^2;
end
%--------------------------------------------------------------------------
% Zero-out boundary rows
D1xo = D1x;
D1yo = D1y;
D0y([1,Ny],:) = 0;
D1y([1,Ny],:) = 0;
D2y([1,Ny],:) = 0;
%--------------------------------------------------------------------------
% Parallel test
if(Nx>1)
    D0x([1,Nx],:) = 0; % Parallel flow test
    D1x([1,Nx],:) = 0; % Parallel flow test
    D2x([1,Nx],:) = 0; % Parallel flow test
end
%--------------------------------------------------------------------------
% Define some useful matrices
Ix = D0x;
Iy = D0y;
blkeye = kron(Ix,Iy);
blkzeros = zeros(Nx*Ny);

% Reshape baseflow
U = reshape(U,Nx*Ny,1);
% dxU = reshape(dxU,Nx*Ny,1);
dyU = reshape(dyU,Nx*Ny,1);
dyyU = reshape(dyyU,Nx*Ny,1);
V = reshape(V,Nx*Ny,1);
dyV = reshape(dyV,Nx*Ny,1);

% % Reshape test
% Ut = U.';
% % dxUt = dxU.';
% dyUt = dyU.';
% dyyUt = dyyU.';
% Vt = V.';
% dyVt = dyV.';
% Ur = zeros(Nx*Ny);
% % dxUr = zeros(Nx*Ny);
% dyUr = zeros(Nx*Ny);
% dyyUr = zeros(Nx*Ny);
% Vr = zeros(Nx*Ny);
% dyVr = zeros(Nx*Ny);
% for i = 1:Nx
%     Ur((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),Ut(i,:));
% %     dxUr((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),dxUt(i,:));
%     dyUr((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),dyUt(i,:));
%     dyyUr((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),dyyUt(i,:));
%     Vr((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),Vt(i,:));
%     dyVr((i-1)*Ny+1:(i-1)*Ny+Ny,:) = kron(ones(Ny,Nx),dyVt(i,:));
% end
% Sponge layer, sigma(x)
sigA = 0.3;
x1 = 49.05; x2 = 850.95;
a1 = sigA/x1^2;
a2 = sigA/(Lx-x2)^2;
sponge = zeros(1,size(xgrid,1));
for i = 1:size(sponge,2)
   if( xgrid(i) <= x1 )
       sponge(i) = a1*(xgrid(i)-x1)^2;
   elseif( xgrid(i) >= x2 )
       sponge(i) = a2*(xgrid(i)-x2)^2;
   end
end
% sigma = kron(diag(sponge),Iy);
sigma = 0.; % For parallel flow case
% figure; plot(xgrid,sponge);

% Laplace operator
% Lap = kron(D2x,Iy)+kron(Ix,D2y)-kz^2*blkeye;
% Rex = Re*delta/delta(end);
Rex = Re;

% K, a common operator
% K = 1/Re*(kron(D2x,Iy) + kron(Ix,D2y) - kz^2*blkeye) + ...
%     - U.*kron(D1x,Iy) - V.*kron(Ix,D1y) - sigma;
K = (kron(1./Rex.*D2x,Iy) + kron(1./Rex.*Ix,D2y) - kz^2*kron(1./Rex.*D0x,D0y)) + ...
    - U.*kron(D1x,Iy) - V.*kron(Ix,D1y) - sigma;
% E matrix (LHS)
E = blkdiag(blkeye,blkeye,blkeye,zeros(Nx*Ny));

% A matrix (RHS): This matrix should be globally stable
A11 = K + dyV.*blkeye;
A12 = -dyU.*blkeye;
A13 = blkzeros;
A14 = -kron(D1x,Iy);
A21 = blkzeros;
% A21 = -dxV.*blkeye; % We have neglected dxV term for this analysis
A22 = K - dyV.*blkeye;
A23 = blkzeros;
A24 = -kron(Ix,D1y);
A31 = blkzeros;
A32 = blkzeros;
A33 = K;
A34 = -zi*kz*blkeye;
A41 = kron(D1xo,eye(Ny));
A42 = kron(eye(Nx),D1yo);
A43 = zi*kz*eye(Nx*Ny);
A44 = blkzeros;
Ao = [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44];

%{
% K, a common operator
K = 1/Re*(kron(D2x,Iy) + kron(Ix,D2y) - kz^2*blkeye) + ...
    - U.*kron(D1x,Iy) - V.*kron(Ix,D1y) - sigma;
% E matrix (LHS)
E = blkdiag(blkeye,blkeye,blkeye,zeros(Nx*Ny));

% A matrix (RHS): This matrix should be globally stable
% A11 = K + dyV.*blkeye;
A11 = K - dxU.*blkeye;
A12 = -dyU.*blkeye;
A13 = blkzeros;
A14 = -kron(D1x,Iy);
A21 = blkzeros;
% A21 = -dxV.*blkeye; % We have neglected dxV term for this analysis
A22 = K - dyV.*blkeye;
A23 = blkzeros;
A24 = -kron(Ix,D1y);
A31 = blkzeros;
A32 = blkzeros;
A33 = K;
A34 = -zi*kz*blkeye;
A41 = kron(D1xo,eye(Ny));
A42 = kron(eye(Nx),D1yo);
A43 = zi*kz*eye(Nx*Ny);
A44 = blkzeros;
Ao = [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44];
%}
% B matrix
Bup = blkdiag(blkzeros,blkzeros,blkzeros);
Bdo = [blkzeros blkzeros blkzeros];
Bo = [Bup; Bdo];
% C matrix
Co = Bo.';

% Impose B.Cs
% Zero-out
%{
% for bi = 1:3
%     for i = 1:Nx
%         E((bi-1)*Nx*Ny+(i-1)*Ny+1,:) = 0;
%         E((bi-1)*Nx*Ny+(i-1)*Ny+Ny,:) = 0;
%         Ao((bi-1)*Nx*Ny+(i-1)*Ny+1,:) = 0;
%         Ao((bi-1)*Nx*Ny+(i-1)*Ny+Ny,:) = 0;
%         Bo((bi-1)*Nx*Ny+(i-1)*Ny+1,:) = 0;
%         Bo((bi-1)*Nx*Ny+(i-1)*Ny+Ny,:) = 0;
%         Co((bi-1)*Nx*Ny+(i-1)*Ny+1,:) = 0;
%         Co((bi-1)*Nx*Ny+(i-1)*Ny+Ny,:) = 0;
%     end
% end
% for bi = 1:3
%     E((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,:) = 0;
%     E((bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny,:) = 0;
%     Ao((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,:) = 0;
%     Ao((bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny,:) = 0;
%     Bo((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,:) = 0;
%     Bo((bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny,:) = 0;
%     Co((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,:) = 0;
%     Co((bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny,:) = 0;
% end
%}
% upper and lower: this task has been done already
for bi = 1:3
    blkind = (bi-1)*Nx*Ny;
    for i = 1:Nx
        subblkind = (i-1)*Ny;
        Ao(blkind+subblkind+1,blkind+subblkind+1:blkind+subblkind+Ny) = [1 zeros(1,Ny-1)];
        Ao(blkind+subblkind+Ny,blkind+subblkind+1:blkind+subblkind+Ny) = [zeros(1,Ny-1) 1];
    end
end
% % left and right (inlet and outlet)
% for bi = 1:3
%     % inlet
%     Ao(blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny,blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny) = eye(Ny);
%     % outlet
%     alp = (xgrid(1)-xgrid(3))/(xgrid(2)-xgrid(3));
%     beta = (xgrid(2)-xgrid(1))/(xgrid(2)-xgrid(3));
%     Ao(blkind+1:blkind+Ny,     blkind+1:blkind+Ny)   = eye(Ny);
%     Ao(blkind+1:blkind+Ny,  blkind+Ny+1:blkind+2*Ny) = -alp*eye(Ny);
%     Ao(blkind+1:blkind+Ny,blkind+2*Ny+1:blkind+3*Ny) = -beta*eye(Ny);
% end
%--------------------------------------------------------------------------
% % Periodic boundary condition test
% if(Nx>1)
%     for bi = 1:3
%         blkind = (bi-1)*Nx*Ny;
%         % inlet
%         Ao(blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny,blkind+Ny+1:blkind+2*Ny) = -eye(Ny);
%         Ao(blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny,blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny) = eye(Ny);
%         % outlet
%         Ao(blkind+1:blkind+Ny,blkind+1:blkind+Ny) = eye(Ny);
%         Ao(blkind+1:blkind+Ny,blkind+(Nx-2)*Ny+1:blkind+(Nx-2)*Ny+Ny) = -eye(Ny);
%     end
% end
%--------------------------------------------------------------------------
% Function output
A = Ao;
B = Bo;
C = Co;
end % End of the function

% state-space form (needed to be further developed with an appropriate
% baseflow)
%{
% E matrix (LHS)
blkeye = kron(Ix,Iy);
% Eo = blkdiag(blkeye,blkeye);
% A matrix (RHS), should be globally stable
A11 = 1/Re*(kron(D4x,Iy)+kron(Ix,D4y)+2*kron(D2x,D2y)-2*kz^2*kron(D2x,Iy)-2*kz^2*kron(Ix,D2y)+kz^4) + ...
    -U.*Lap*kron(D1x,Iy) + ...
    -V.*Lap*kron(Ix,D1y) + ...
    -dV.*Lap + ...
    -2*dU.*kron(D2x,Iy) + ...
    -d2V.*kron(Ix,D1y) + ...
    +d2U.*kron(D1x,Iy) + ...
    -d3V + ...
    -2*(dxyU.*kron(D1x,Iy) + dxU.*kron(Dx1,Dy1)/(kron(D2x,Iy) - kz^2)*kron(Dx1,Dy1));
A11 = Lap\A11 - sigma;
A12 = 2*zi*kz*(dxyU.*kron(D1x,Iy) + dxU.*kron(D1x,D1y))/(kron(D2x,Iy) - kz^2);
A21 = -zi*kz*dyU*blkeye;
A22 = 1/Re*Lap - U.*kron(D1x,Iy) + ...
    -V.*kron(Ix,D1y) + ...
    -dxU.*blkeye;
A22 = Lap\A22 - sigma;
Ao   = [A11 A12; A21 A22];

% B
B11 = -Lap\(f.*kron(D1x,D1y) + df.*kron(D1x,Iy));
B12 = Lap\(f.*kron(D2x,Iy) - kz^2*f*blkeye);
B13 = -Lap\zi*kz*(df*blkeye + f.*kron(Ix,D1y));
B21 = zi*kz*f*blkeye;
B22 = zeros(Nx*Ny);
B23 = -f.*kron(D1x,Iy);
Bo = [B11 B12 B13; B21 B22 B23];

% C
Cio = (kron(D2x,Iy) - kz^2); % inversed operator
C11 = -Cio\kron(D1x,D1y);
C12 = Cio\zi*kz*blkeye;
C21 = blkeye;
C22 = zeros(Nx*Ny);
C31 = -Cio\zi*kz*kron(Ix,D1y);
C32 = -Cio\kron(D1x,Iy);
Co = [C11 C12; C21 C22; C31 C32];

%         C = blkdiag(D01,D01,D01)\C;
% Function output
A = Ao;
B = Bo;
C = Co;
% end
%}