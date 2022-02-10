function w = arbSpaceWeight(M,N,x0,xgrid)
%--------------------------------------------------------------------------
% Arbitary grid spacing weight for PS derivative matrix including M = 0.
% Input -------------------------------------------------------------------
% M: order of the derivatives
% N: number of grids
% x0: evaluating location
% xgrid: collocation grids
% Output ------------------------------------------------------------------
% w: M x N weight matrix for xgrid.
%--------------------------------------------------------------------------
w = zeros(N,N,M+1); % n, nu, m
w(1,1,1) = 1; % d0,00
c1 = 1;
for n = 2:N
    c2 = 1;
    for nu = 1:n-1
        c3 = xgrid(n)-xgrid(nu);
        c2 = c2*c3;
%         if(n<=M)
%             dnn_1nu = 0;
%         end
        for m = 0:min(n-1,M)
            if(m==0)
                w(n,nu,m+1) = (xgrid(n)-x0)*w(n-1,nu,m+1)/c3;
            else
                w(n,nu,m+1) = ((xgrid(n)-x0)*w(n-1,nu,m+1)-(m)*w(n-1,nu,m))/c3;
            end
        end % End of m for loop
    end % End of nu for loop 
    for m = 0:min(n-1,M)
        if(m==0)
            w(n,n,m+1) = c1/c2*(-(xgrid(n-1)-x0)*w(n-1,n-1,m+1));
        else
            w(n,n,m+1) = c1/c2*((m)*w(n-1,n-1,m)-(xgrid(n-1)-x0)*w(n-1,n-1,m+1));
        end
    end % End of m for loop
    c1 = c2;
end % End of n for loop
end % End of the function