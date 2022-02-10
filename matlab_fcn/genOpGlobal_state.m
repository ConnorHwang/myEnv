function [A,B,C,E] = genOpGlobal_state(Lx,Nx,Ly,Ny)
global kz kx zi Re U dyU dyyU V dyV
% generateOperators
[~,D0x,D1x,D2x,D3x,D4x] = genChebGlobal(Nx,Lx);
[~,D0y,D1y,D2y,~,D4y] = genChebGlobal(Ny,Ly);
%--------------------------------------------------------------------------
% Parallel test
if( Nx == 1 )
    D1x = zi*kx;
    D2x = -kx^2;
    D3x = -zi*kx^3;
    D4x =  kx^4;
end
%--------------------------------------------------------------------------
% Zero-out boundary rows
% D1xo = D1x;
% D1yo = D1y;
D0y([1,Ny],:) = 0;
D0y([2,Ny-1],:) = 0;
% D1y([1,Ny],:) = 0;
D2y([1,Ny],:) = 0;
D2y([2,Ny-1],:) = 0;
% D3y([1,Ny],:) = 0;
D4y([1,Ny],:) = 0;
D4y([2,Ny-1],:) = 0;
%--------------------------------------------------------------------------
% Parallel test
if( Nx > 1 )
    D0x([1,Nx],:) = 0; % Parallel flow test
    D1x([1,Nx],:) = 0; % Parallel flow test
    D2x([1,Nx],:) = 0; % Parallel flow test
    D3x([1,Nx],:) = 0; % Parallel flow test
    D4x([1,Nx],:) = 0; % Parallel flow test
end
%--------------------------------------------------------------------------
% Define some useful matrices
Ix = D0x;
Iy = D0y;

% Reshape baseflow
U = reshape(U,Nx*Ny,1);
% dxU = reshape(dxU,Nx*Ny,1);
dyU = reshape(dyU,Nx*Ny,1);
dyyU = reshape(dyyU,Nx*Ny,1);
V = reshape(V,Nx*Ny,1);
dyV = reshape(dyV,Nx*Ny,1);

% figure; plot(xgrid,sponge);

% Laplace operator
% Lap = kron(D2x,Iy)+kron(Ix,D2y)-kz^2*blkeye;
% Rex = Re*delta/delta(end);
Rex = Re;

% Explicit boundary condition with Weideman's operators
        Lap = kron(D2x,Iy) + kron(Ix,D2y) - kz^2*kron(Ix,Iy);
        % M
        M11 = Lap;
        M22 = kron(Ix,Iy);
        M   = blkdiag(M11,M22);
        % A
        A11 = 1./Rex.*(kron(D4x,Iy) + kron(Ix,D4y) + kz^4*kron(Ix,Iy) +...
              2*(kron(D2x,D2y) - kz^2*kron(D2x,Iy) - kz^2*kron(Ix,D2y))) + ...
            - U.*(kron(D3x,Iy) + kron(D1x,D2y) - kz^2*kron(D1x,Iy)) + ...
              dyyU.*kron(D1x,Iy);
        A12 = zeros(Nx*Ny);
        A21 =-zi*kz*dyU.*kron(Ix,Iy);
        A22 = 1./Rex.*Lap - U.*kron(D1x,Iy);
        A   = [A11 A12; A21 A22];

        % Divide by an LHS operator
%         Ao = M\A;
%         B = M\B;
%         C = blkdiag(D01,D01,D01)\C;
        Eo = zeros(Nx*Ny*2);        
        Bo = zeros(Nx*Ny*2);
        Co = zeros(Nx*Ny*2);

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
for bi = 1:2
    blkind = (bi-1)*Nx*Ny;
    for i = 1:Nx
        subblkind = (i-1)*Ny;
        M(blkind+subblkind+1,blkind+subblkind+1:blkind+subblkind+Ny) = [1 zeros(1,Ny-1)];
        M(blkind+subblkind+Ny,blkind+subblkind+1:blkind+subblkind+Ny) = [zeros(1,Ny-1) 1];
        M(blkind+subblkind+2,blkind+subblkind+1:blkind+subblkind+Ny) = D1y(1,:);
        M(blkind+subblkind+Ny-1,blkind+subblkind+1:blkind+subblkind+Ny) = D1y(Ny,:);
    end
end
% A(2,:) = 0;
% A(Ny-1,:) = 0;
% A(Ny+2,:) = 0;
% A(2*Ny-1,:) = 0;
% % left and right (inlet and outlet)
% for bi = 1:3
%     % inlet
%     Ao((bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny,(bi-1)*Nx*Ny+(Nx-1)*Ny+1:(bi-1)*Nx*Ny+(Nx-1)*Ny+Ny) = eye(Ny);
%     % outlet
%     alp = (xgrid(1)-xgrid(3))/(xgrid(2)-xgrid(3));
%     beta = (xgrid(2)-xgrid(1))/(xgrid(2)-xgrid(3));
%     Ao((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,     (bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny)   = eye(Ny);
%     Ao((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,  (bi-1)*Nx*Ny+Ny+1:(bi-1)*Nx*Ny+2*Ny) = -alp*eye(Ny);
%     Ao((bi-1)*Nx*Ny+1:(bi-1)*Nx*Ny+Ny,(bi-1)*Nx*Ny+2*Ny+1:(bi-1)*Nx*Ny+3*Ny) = -beta*eye(Ny);
% end
%--------------------------------------------------------------------------
% Periodic boundary condition test
if( Nx > 1 )
    disp('Nx > 1, Periodic boundary condition')
    for bi = 1:2
        blkind = (bi-1)*Nx*Ny;
        % inlet
        M(blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny,blkind+Ny+1:blkind+2*Ny) = -eye(Ny);
        M(blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny,blkind+(Nx-1)*Ny+1:blkind+(Nx-1)*Ny+Ny) = eye(Ny);
        % outlet
        M(blkind+1:blkind+Ny,blkind+1:blkind+Ny) = eye(Ny);
        M(blkind+1:blkind+Ny,blkind+(Nx-2)*Ny+1:blkind+(Nx-2)*Ny+Ny) = -eye(Ny);
    end
end
%--------------------------------------------------------------------------
% Function output
A = M\A;
B = Bo;
C = Co;
E = Eo;
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