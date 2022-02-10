function c = intMat(N,Lmin,Lmax,coord)
% Integration weighting matrix
% input: A state vector (either in cylindrical or Cartisian coordinate)
% output: \int_domain 1/2*(u_hat_i^2) dy.
% Check if you have to do Cholesky decomposition.
if(strcmp(coord,'cart'))
    num = round(abs(N));
    c = zeros(num,num);
    L = Lmax-Lmin;
    for i = (0:num-1)
        for j = (0:num-1)
            if rem(i+j,2) == 0
                p = (1/(1-(i+j)^2)+1/(1-(i-j)^2))*(L/2); % 1/2 due to the linear mapping
                c(i+1,j+1)=p;
            else
                c(i+1,j+1)=0;
            end
        end
    end
elseif(strcmp(coord,'cyl')) % Cylindrical coordinate
    % Returns L2 norm of an state vector in Chebyshev space (Gauss-Lobatto)
    num = round(abs(N));
    c = zeros(num,num);
    a = (Lmax - Lmin)/2;
    b = (Lmin + Lmax)/2;
    
    for i = (0:num-1)
        for j = (0:num-1)
            if rem(i+j,2) == 0
                p = a*b*(-((1+(-1)^(i+j))*(-1+i^2+j^2))/(i^4+(-1+j^2)^2-2*i^2*(1+j^2)));
                c(i+1,j+1)=p;
            else
                p = a*a*((-1+(-1)^(i+j))*(-4+i^2+j^2))/(i^4+(-4+j^2)^2-2*i^2*(4+j^2));
                c(i+1,j+1)=p;
            end
        end
    end
    % c_ = chol(c);
end % End of if

end % End of the function