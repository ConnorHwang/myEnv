function DM = genFPOp(x,M,N)
% Generate derivative operators using arbitrary spacing algorithm
% (arbSpaceWeight.m) in order to use fictitious point (FP) method
% Input -------------------------------------------------------------------
% a  : Discretized grids
% xfp: Fictitious point
% M  : The highest number of derivatives to print out
% N  : 
% Output ------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Container
DM = zeros(N,N,M+1);
% Fictitious point method
% Define new collocation points
% newx = [x0; rp1];
% N = size(newx,1);
% D0t = zeros(N,N);
% D1t = zeros(N,N);
% D2t = zeros(N,N);
% D3t = zeros(N,N);
%         t = arbSpaceWeight(M,N,x0,newx);
% Build system
for i = 1:N
    xeval = x(i);
    t = arbSpaceWeight(M,N,xeval,x);
    for j = 1:M+1
        DM(i,:,j) = t(N,:,j);
    end
end
% 
% % Create boundary condition
% BCt = arbSpaceWeight(M,N,rp1(end),newx);
% BC = BCt(end,:,2);
% % Subtract BC and make zero column for x_fp weight
% for i = 1:N-1
%     D3t(i+1,:) = D3t(i+1,:) - D3t(i+1,1)/BC(1,1)*BC;
% end
% % Remove x_fp row
% D3t(1,:) = [];
% D3t(:,1) = [];
% % Impose other boundary conditions
% D3t([1,N1],:) = [];
% D3t(:,[1,N1]) = [];
% Nm = size(D3t,1);
% A1 = D3t;
% E1 = eye(Nm);
% %         disp('finish')
end % End of the function