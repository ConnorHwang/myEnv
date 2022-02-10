function cr = recepCoeff(phi,B,rp1,del99,option)
global N1 L1
switch(option)
    case 1 % 2D parallel flow
%         Iwblck = sqrt(ChebQuad(N1+2,L1)); Iwblck = diag(Iwblck(2:N1+1));
%         Dgs = diag(rp1<del99)*Iwblck;
        Dgs = diag(rp1<del99);
        Dg = blkdiag(Dgs,Dgs,Dgs);
%         Wg = Dg*W;
%         Wg = W;
%         Wg = eye(3*N1);
%         cr_ = trace(Dg*phi*Dg')/trace(B*Dg*B');
%         cr_ = trace(phi)/trace(B*Dg*B');
        cr_ = trace(phi)/trace(B);
%         log10(cr_)
        if( abs(imag(cr_)) > 10^-6)
            fprintf('log10(Cr) = %f (Cr = %f + %fi)\n', log10(cr_), real(cr_), imag(cr_));
            fprintf('Imaginary part of Cr is larger than the tolerance!\n');
        end
        cr = real(cr_);
%         log10(cr)
%         fprintf('log10(Cr) = %f (Cr = %f + %fi)\n', log10(cr), real(cr), imag(cr));
%         if( ~isreal(cr) )
%         if( cr<0 )
% %             fprintf('Receptive coefficient is not real!\n');
%             fprintf('log10(Cr) = %f (%f + %fi)\n', log10(cr), real(cr), imag(cr));
%         end
        % test
%         a = real(trace(Dg*phi*Dg'));
%         log10(a/2.0415)
%         test = f(rp1<5);
    case 2 % Global
end % End of 'switch'
end