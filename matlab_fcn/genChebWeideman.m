function [x,D0,D1,D2,D3,D4] = genChebWeideman(N,L)
% using cheb4c.m and chebdif.m. Note that boundary conditions are already
% imbedded and check the matrix size.
D0 = eye(N);
[~,DM] = chebdif(N+2,3);
D1_ = DM(2:N+1,2:N+1,1);
D2_ = DM(2:N+1,2:N+1,2);
D3_ = DM(2:N+1,2:N+1,3);
[x,D4_] = cheb4c(N+2);
% re-scaling Chebyshev points
x = (x+1)/2*L;
% re-scaling Chebyshev operators
dy_cheb_dy_phy = 2/L;
D1 = D1_*(dy_cheb_dy_phy);
D2 = D2_*(dy_cheb_dy_phy)^2;
D3 = D3_*(dy_cheb_dy_phy)^3;
D4 = D4_*(dy_cheb_dy_phy)^4;

end % End of function