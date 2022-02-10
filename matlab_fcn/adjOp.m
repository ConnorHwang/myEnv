function B = adjOp(A,N,Lmin,Lmax)
% returns adjoint of an operator A
Astar = conj(A);
% <h2,Q2*A*h1> = <Astar*h2,Q1*h1>, the point is weighting should go into
% the middle of the two terms that are being producted. (e.g., h1*Q*h2)
% assuming integration weights in H1 and H2 for h1 and h2. And H1 and H2
% are the same.
Q = L2Cheb(N,Lmin,Lmax);
B = Q\Astar*Q;
% Boundary conditions for B?
end