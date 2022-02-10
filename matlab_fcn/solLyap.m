function c = solLyap(a,b)
c = lyap(a,b);
diag_c = diag(c);
if( sum(diag_c<0) )
    fprintf('Lyapunov solution is not positive definite!\n');
end
end % End of function