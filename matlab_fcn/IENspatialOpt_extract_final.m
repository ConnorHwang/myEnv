function [M_growth, initial, final, Kappa, Q0origin, cols] = ...
    IENspatialOpt_extract_final(We, n, r, N1, L1, N2, L2, Eigenvectors, Eigenvalues, K, t, normtype)
%--------------------------------------------------------------------------
% Brief : Primitive variable code, non-modal growth with seminorm. 
% Input : 
% 'normType' - (0) kinetic energy of the field, (1): interface.
% 'K' - # of modes in consider.
% 'Eigenvectors' - Entire eigenvectors.
% 'Eigenvalues' - Entire eigenvalues.
% Output: Transient growth rate. Initial(right singular vector)
%         and final(left singular vector).
%--------------------------------------------------------------------------
zi = sqrt(-1);
% Define vectors following main code's notation.
V2vec = 1:N2;
P2vec = 3*N2+1:4*N2;
V2_vec = 4*N2+1:5*N2;
W2_vec = 5*N2+1:6*N2;
U2_vec = 6*N2+1:7*N2;
Interface = 7*N2+1;
Interface_ = 7*N2+1+1;
ublk = Interface_; % Size of the upper matrix block
V1vec = ublk+1:ublk+N1;
P1vec = ublk+3*N1+1:ublk+4*N1;
V1_vec = ublk+4*N1+1:ublk+5*N1;
W1_vec = ublk+5*N1+1:ublk+6*N1;
U1_vec = ublk+6*N1+1:ublk+7*N1;

% Use as many as modes if the number of eigenvalue is less than the
% specified value 'K'.
if( size(Eigenvectors,2) < K )
%     fprintf('Warning: K = %g is used instead of specified value K = %g.\n', size(Eigenvectors,2), K);
    K = size(Eigenvectors,2);
end
% Selecting the number of Eigenmodes user want to consider
ishift = 1;
% This while loop excludes the unstable modes.
% while imag(Eigenvalues(ishift))<=0, ishift=ishift+1; end % note for temporal code >=.
cols = (ishift:ishift + K - 1);
% Building Q0 matrix
Q0 = Eigenvectors(:, cols);
% To recover the ansatz of the eigenvector, multiply zi for radial comp.
Q0([V2vec V1vec],:) = Q0([V2vec V1vec],:) * zi;
% Save Eigenvectors before eliminating the pressure terms
Q0origin = Q0;
%--------------------------------------------------------------------------
% Q1root   = Q0;
% By commenting out this section, we are not using weighting
% Multiply weight for interface distortion, B = abs(alp^2+n^2-1)/We*2.
% B_ = sqrt((Eigenvalues(cols).^2+n^2)/We*2).';
% B_ = sqrt((-1+Eigenvalues(cols).^2+n^2)/We).'; %test
B_ = sqrt((-1+Eigenvalues(cols).^2+n^2)/We).'; %test
Q0(Interface,:) = B_ .* Q0(Interface,:);
% Q1root([P2vec V2_vec W2_vec U2_vec Interface_ P1vec V1_vec W1_vec U1_vec],:) = [];
%--------------------------------------------------------------------------
% Get rid of Pressure temrs. For spatial code, get rid of extra dummy
% terms.
Q0([P2vec V2_vec W2_vec U2_vec Interface_ P1vec V1_vec W1_vec U1_vec],:) = [];
% Build peudo fundamental operator, Mjet function give the matrix of which
% multiplication to a state vector 'q0' will give sqrt(2 norm).
[~, M_int_1] = MJet(N1,0,L1);
[~, M_int_2] = MJet(N2,L1,L1+L2);
M_int_2 = M_int_2 * sqrt(r);
M_interface = 1;
M_int_1 = M_int_1 * sqrt(1/2); % test
M_int_2 = M_int_2 * sqrt(1/2); % test
% M_interface = sqrt(B);
% Constraints to the all space.
CONSTRAINT = blkdiag(M_int_2,M_int_2,M_int_2,M_interface,...
                     M_int_1,M_int_1,M_int_1);
% Customized energy Matrix
IEN = blkdiag(zeros(N2), zeros(N2), zeros(N2),M_interface,...
                 zeros(N1), zeros(N1), zeros(N1));
% V-norm
VNORM = blkdiag(M_int_2, zeros(N2), zeros(N2), 0,...
                 M_int_1, zeros(N1), zeros(N1));
% Connor's method,
Fc = CONSTRAINT * Q0; % which has a full column rank.
Fcp = pinv(Fc);
if ( normtype == 0 ) % normtype = 0 for the regular norm.
    qb = zi*Fc*diag(Eigenvalues(cols))*Fcp; % For temporal, follow -i*omega*t;
    svd_ = expm(t*qb);
    [~,S,V] = svd(svd_);
    M_growth = S(1,1)^2;
    % Evolve Q0_1
    Qfinal = Q0origin .* exp(zi*(Eigenvalues(cols)).'* t);
    Kappa =  Fcp * V(:,1);
    initial = Q0origin * Kappa;
    final = Qfinal * Kappa;
elseif ( normtype == 1 ) % normtype = 1 for the semi norm.
    if (t ~= 0)
        Q1 = Q0;
%         Q1 = Q1root; % CHECK
        for i = cols
           Q1(:,i) = exp(zi*t*Eigenvalues(i))*Q1(:,i); 
        end
%         Q1(Interface,:) = B_ .* Q1(Interface,:); % CHECK
        svd_ = IEN * Q1 * Fcp;
        % M_growth = norm(expm(svd_))^2;
        [~,S,V] = svd(svd_);
        M_growth = S(1,1)^2;
        Qfinal = Q0origin .* exp(zi*(Eigenvalues(cols)).'* t);
        Kappa =  Fcp * V(:,1);
        initial = Q0origin * Kappa;
        final = Qfinal * Kappa;
    else 
%         fprintf('t = 0!\n');
        svd_ = (IEN*Q0)*Fcp;
        [~,S,V] = svd(svd_);
        M_growth = S(1,1)^2;
        % Evolve Q0_1
        Qfinal = Q0origin;
        Kappa =  Fcp * V(:,1);
        initial = Q0origin * Kappa;
        final = Qfinal * Kappa;
    end
else % normtype = 2 for v-norm
    if (t ~= 0)
        Q1 = Q0;
        for i = cols
           Q1(:,i) = exp(zi*t*Eigenvalues(i))*Q1(:,i); 
        end
        svd_ = VNORM * Q1 * Fcp;
        % M_growth = norm(expm(svd_))^2;
        [~,S,V] = svd(svd_);
        M_growth = S(1,1)^2;
        Qfinal = Q0origin .* exp(zi*(Eigenvalues(cols)).'* t);
        Kappa =  Fcp * V(:,1);
        initial = Q0origin * Kappa;
        final = Qfinal * Kappa;
    else 
%         fprintf('t = 0!\n');
        svd_ = (VNORM*Q0)*Fcp;
        [~,S,V] = svd(svd_);
        M_growth = S(1,1)^2;
        % Evolve Q0_1
        Qfinal = Q0origin;
        Kappa =  Fcp * V(:,1);
        initial = Q0origin * Kappa;
        final = Qfinal * Kappa;
    end
end % End of semi-norm.
end