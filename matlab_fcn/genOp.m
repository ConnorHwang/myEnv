function [A,B,C] = genOp(N1,option)
global D01 D11 D21 D41 kz kx zi Re u1 du1 d2u1 f df
% generateOperators
switch(option)
    case 1
        % Option 1, normal velocity and vormal vorticity, 
        % Explicit boundary condition with classical D matrices
        k2  = kx^2+kz^2;
        Lap = D21-k2*D01;
        % M
        M11 = Lap;
        M22 = D01;
        M   = blkdiag(M11,M22);
        % A
        A11 = 1/Re*(D41-2*k2*D21+k2^2*D01) - zi*kx*u1.*Lap + zi*kx*d2u1.*D01;
        A12 = zeros(N1);
        A21 =-zi*kz*du1.*D01;
        A22 = 1/Re*Lap - zi*kx*u1.*D01;
        A   = [A11 A12; A21 A22];
        
        % B.C.
        M = imposeBC(M,N1,1,1, 1,0, 1);
        M = imposeBC(M,N1,1,1,-1,0,-1);
        M = imposeBC(M,N1,1,1, 2,1, 1);
        M = imposeBC(M,N1,1,1,-2,1,-1);
        M = imposeBC(M,N1,2,2, 1,0, 1);
        M = imposeBC(M,N1,2,2,-1,0,-1);

%         % B.C. zeros
%         A = imposeBC(A,N1,1,1,[1,-1],-1, 0);
%         A = imposeBC(A,N1,2,2,[1,-1],-1, 0);
        % B.C. er * D
        A = imposeBC(A,N1,1,1, 1,2, 1);
        A = imposeBC(A,N1,1,1,-1,2,-1);
        A = imposeBC(A,N1,1,1, 2,3, 1);
        A = imposeBC(A,N1,1,1,-2,3,-1);
        A = imposeBC(A,N1,2,2, 1,2, 1);
        A = imposeBC(A,N1,2,2,-1,2,-1);
%         % B.C. D
%         A = imposeBC(A,N1,1,1, 1,0, 1);
%         A = imposeBC(A,N1,1,1,-1,0,-1);
%         A = imposeBC(A,N1,1,1, 2,1, 1);
%         A = imposeBC(A,N1,1,1,-2,1,-1);
%         A = imposeBC(A,N1,2,2, 1,0, 1);
%         A = imposeBC(A,N1,2,2,-1,0,-1);
        % B
        B11 = -(zi*kx*f.*D11 + zi*kx*df.*D01); %
        B12 = -(kx^2*f.*D01 + kz^2*f.*D01); %
        B13 = -zi*kz*(df.*D01 + f.*D11); %
        B21 = zi*kz*f.*D01; % sign error, should be positive?
        B22 =  zeros(N1);
        B23 = -zi*kx*f.*D01;
        B = [B11 B12 B13; B21 B22 B23];
        % B.C.
        B = imposeBC(B,N1,[1,2],[1,2,3],[1,-1],-1,0);
        B = imposeBC(B,N1,1    ,[1,2,3],[2,-2],-1,0);
        % C
        C11 =  k2\zi*kx*D11;
        C12 = -k2\zi*kz*D01;
        C21 = eye(N1)*D01;
        C22 = zeros(N1);
        C31 =  k2\zi*kz*D11;
        C32 =  k2\zi*kx*D01;
        C = [C11 C12; C21 C22; C31 C32];
        % B.C.
%         C = imposeBC(C,N1,[1,2,3],[1,2], 1, 0, 1);
%         C = imposeBC(C,N1,[1,2,3],[1,2],-1, 0,-1);
        C = imposeBC(C,N1,[1,2,3],1    , 2, 1, 1);
        C = imposeBC(C,N1,[1,2,3],1    ,-2, 1,-1);
        C = imposeBC(C,N1,[1,2,3],[1,2],[1,-1],-1, 0);
        C = imposeBC(C,N1,[1,2,3],1    ,[1,-1],-1, 0);
        % Divide by an LHS operator
        A = M\A;
        B = M\B;
%         C = blkdiag(D01,D01,D01)\C;
    case 2
        % Using Weideman's operator
        k2  = kx^2+kz^2;
        Lap = D21-k2*D01;
        % M
        M11 = Lap;
        M22 = D01;
        M   = blkdiag(M11,M22);
        % A
        A11 = 1/Re*(D41-2*k2*D21+k2^2*D01)-zi*kx*u1.*Lap+zi*kx*d2u1.*D01;
        A12 = zeros(N1);
        A21 = -zi*kz*du1.*D01;
        A22 = 1/Re*Lap-zi*kx*u1.*D01;
        A   = [A11 A12; A21 A22];
        
        % B
        B11 = -(zi*kx*f.*D11 + zi*kx*df.*D01); %
        B12 = -(kx^2*f.*D01 + kz^2*f.*D01); %
        B13 = -zi*kz*(df.*D01 + f.*D11); %
        B21 =  zi*kz*f.*D01; % sign error, should be positive?
        B22 =  zeros(N1);
        B23 = -zi*kx*f.*D01;
        B   = [B11 B12 B13; B21 B22 B23];
        
        % C
        C11 =  k2\zi*kx*D11;
        C12 = -k2\zi*kz*D01;
        C21 =  eye(N1)*D01;
        C22 =  zeros(N1);
        C31 =  k2\zi*kz*D11;
        C32 =  k2\zi*kx*D01;
        C   = [C11 C12; C21 C22; C31 C32];
        
        % Divide by an LHS operator
        A = M\A;
        B = M\B;
    case 3
        % Explicit boundary condition with Weideman's operators
        k2  = kx^2+kz^2;
        Lap = D21-k2*D01;
        % M
        M11 = Lap;
        M22 = D01;
        M   = blkdiag(M11,M22);
        % A
        A11 = 1/Re*(D41-2*k2*D21+k2^2*D01) - zi*kx*u1.*Lap + zi*kx*d2u1.*D01;
        A12 = zeros(N1);
        A21 =-zi*kz*du1.*D01;
        A22 = 1/Re*Lap - zi*kx*u1.*D01;
        A   = [A11 A12; A21 A22];
        
        % B.C.
        M = imposeBC(M,N1,1,1, 1,0, 1);
        M = imposeBC(M,N1,1,1,-1,0,-1);
        M = imposeBC(M,N1,1,1, 2,1, 1);
        M = imposeBC(M,N1,1,1,-2,1,-1);
        M = imposeBC(M,N1,2,2, 1,0, 1);
        M = imposeBC(M,N1,2,2,-1,0,-1);

%         % B.C. zeros
%         A = imposeBC(A,N1,1,1,[1,-1],-1, 0);
%         A = imposeBC(A,N1,2,2,[1,-1],-1, 0);
        % B.C. er * D
        A = imposeBC(A,N1,1,1, 1,2, 1);
        A = imposeBC(A,N1,1,1,-1,2,-1);
        A = imposeBC(A,N1,1,1, 2,3, 1);
        A = imposeBC(A,N1,1,1,-2,3,-1);
        A = imposeBC(A,N1,2,2, 1,2, 1);
        A = imposeBC(A,N1,2,2,-1,2,-1);
%         % B.C. D
%         A = imposeBC(A,N1,1,1, 1,0, 1);
%         A = imposeBC(A,N1,1,1,-1,0,-1);
%         A = imposeBC(A,N1,1,1, 2,1, 1);
%         A = imposeBC(A,N1,1,1,-2,1,-1);
%         A = imposeBC(A,N1,2,2, 1,0, 1);
%         A = imposeBC(A,N1,2,2,-1,0,-1);
        % B
        B11 = -(zi*kx*f.*D11 + zi*kx*df.*D01); %
        B12 = -(kx^2*f.*D01 + kz^2*f.*D01); %
        B13 = -zi*kz*(df.*D01 + f.*D11); %
        B21 = zi*kz*f.*D01; % sign error, should be positive?
        B22 =  zeros(N1);
        B23 = -zi*kx*f.*D01;
        B = [B11 B12 B13; B21 B22 B23];
        % B.C.
        B = imposeBC(B,N1,[1,2],[1,2,3],[1,-1],-1,0);
        B = imposeBC(B,N1,1    ,[1,2,3],[2,-2],-1,0);
        % C
        C11 =  k2\zi*kx*D11;
        C12 = -k2\zi*kz*D01;
        C21 = eye(N1)*D01;
        C22 = zeros(N1);
        C31 =  k2\zi*kz*D11;
        C32 =  k2\zi*kx*D01;
        C = [C11 C12; C21 C22; C31 C32];
        % B.C.
%         C = imposeBC(C,N1,[1,2,3],[1,2], 1, 0, 1);
%         C = imposeBC(C,N1,[1,2,3],[1,2],-1, 0,-1);
        C = imposeBC(C,N1,[1,2,3],1    , 2, 1, 1);
        C = imposeBC(C,N1,[1,2,3],1    ,-2, 1,-1);
        C = imposeBC(C,N1,[1,2,3],[1,2],[1,-1],-1, 0);
        C = imposeBC(C,N1,[1,2,3],1    ,[1,-1],-1, 0);
        % Divide by an LHS operator
        A = M\A;
        B = M\B;
%         C = blkdiag(D01,D01,D01)\C;
end % End of 'switch'
end