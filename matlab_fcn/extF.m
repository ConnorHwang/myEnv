function [ef,W] = extF(B,N1,option)
% global L1
switch(option)
    case 1
%         Q  = energyMat(N1,L1);
%         Iwblck = chol(intMat(N1,0,L1,'cart'));
%         test = Iwblck*Iwblck';
%         W = blkdiag(test,test,test);
        W = eye(3*N1); % 1, white-in-space noise; 2, HIT
        ef = B*W*B'; % Note that the function 'lyap' solves A*X + X*A' + Q = 0
%         diag_ef = diag(ef);
        % Check positive definite
%         [~,pd] = chol(ef);
%         if( pd > 0 )
%            fprintf('Forcing term is not positive definite!\n'); 
%         end
    case 2
end % End of 'switch'
end