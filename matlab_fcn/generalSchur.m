function [R,L] = generalSchur(A,B,C,D,E,F,option)
switch(option)
    case 1
        % AR-LB=C, DR-LE=F
        % Solve generalized Sylvester equation using Kronecker product notation
        p = size(A,1);
        q = size(B,1);
        % lhs = [kron(eye(p),A) kron(-B.',eye(q)); kron(eye(p),D) kron(-E.',eye(q))];
        lhs = [kron(eye(q),A) kron(-B.',eye(p)); kron(eye(q),D) kron(-E.',eye(p))];
        rhs = [reshape(C,[p*q,1]); reshape(F,[p*q,1])];
        tic;
        x   = lhs\rhs;
        time = toc;
        X   = reshape(x,[p,2*q]);
        R   = X(:,1:q);
        L   = X(:,q+1:2*q);
        
        checkC = A*R - L*B;
        checkF = D*R - L*E;
        
        cerror = sum(sum(checkC - C));
        ferror = sum(sum(checkF - F));
        fprintf('Time for solving (Kronecker): %f [sec]\n', time);
        tol = 1e-6;
        if((abs(cerror) > tol)||(abs(ferror) > tol))
            fprintf('C error: %f\nF error: %f\n', abs(cerror), abs(ferror));
        end
        
    case 2
        % Solves the genralized Schur equation following the algorithm introduced
        % in Kagtrom & Westin(1989), IEEE
        % AR-LB=C; DR-LE=F;
        % (A,D): p by p, (B,E): q by q matrix, (R,L,C,F): p by q matrix
        tic; % time measure
        p = size(A,1);
        q = size(B,1);
        R1 = zeros(p,q);
        L1 = zeros(p,q);
        % STEP 1: Transform (A,D) and (B,E) into simpler form
        % A1 is upper Hessenberg
        % D1 is upper triangular
        % B1 is upper quasi-triangular
        % E1 is upper triangular
        % upper hessenberg matrix AA, upper triangular matrix BB
        % unitary matrices Q and Z such that, Q*A*Z = AA and Q*B*Z = BB
        [A1,D1,Qad,Zad] = hess(A,D);
        % [A1,D1,Qad,Zad] = qz(A,D);
        [B1,E1,Qbe,Zbe] = qz(B,E);
        % [A1,D1,Qad,Zad] = ordqz(A1,D1,Qad,Zad,'rhp');
        % [B1,E1,Qbe,Zbe] = ordqz(B1,E1,Qbe,Zbe,'rhp');
        
        % STEP 2: Modify the right-hand sides (C,F): C1 = P^T*C*V; F1 = P^T*F*V
        C1 = Qad*C*Zbe;
        F1 = Qad*F*Zbe;
        
        % STEP 3-(1): Solve the transformed system for L1 and R1
        % GS-Algorithm:
        for j = 1:q
            for i = p:-1:1
                Abuf = [A1(i,i) -B1(j,j); D1(i,i) -E1(j,j)];
                Bbuf = [C1(i,j); F1(i,j)];
                Xbuf = Abuf\Bbuf;
                R1(i,j) = Xbuf(1);
                L1(i,j) = Xbuf(2);
                for k = 1:(i-1)
                    C1(k,j) = C1(k,j)-A1(k,i)*R1(i,j);
                    F1(k,j) = F1(k,j)-D1(k,i)*R1(i,j);
                end
                for k = (j+1):q
                    C1(i,k) = C1(i,k)+L1(i,j)*B1(j,k);
                    F1(i,k) = F1(i,k)+L1(i,j)*E1(j,k);
                end
            end
        end % End of 'j' for loop
        
        % STEP 4: Transform the solution back to the original system. L = P^T*L1*U,
        % R = Q*R1*V^T
        L = Qad*L1*Qbe.';
        R = Zad*R1*Zbe.';
        
        
        
        time = toc;
        % STEP 3-(2): Check your soln.
        checkC = A*R - L*B;
        checkF = D*R - L*E;
        
%         cerror = sum(sum(abs((checkC - C)./C)));
%         ferror = sum(sum(abs((checkF - F)./F)));
        cerror = sum(sum(abs((checkC - C))));
        ferror = sum(sum(abs((checkF - F))));
        fprintf('Time for solving (GS-algorithm): %f [sec]\n', time);
        tol = 1e-6;
        if((abs(cerror) > tol)||(abs(ferror) > tol))
            fprintf('GS alogrithm\n');
            fprintf('C error: %f\nF error: %f\n', cerror, ferror);
        end
    case 3
        % Solves the genralized Schur equation following the algorithm introduced
        % in A. Shahzad et al. Automatica 47 (2011) 244--248.
        % AR-LB=C; DR-LE=F;
        % Conditions for the input matrix
        % A = pxp, upper-triangular matrix
        % B = qxq, all the diagonal eliments must be zero, upper-triangular matrix
        % A and D matrices become upper-triangular by qz algorithm
        % q >= 1. This requirement automatically satisfied if E is a singular matrix
        
        % Check the requirements
        detA = det(A);
        detE = det(E);
        diagB = sum(abs(diag(B)));
        fprintf('det(A) = %e, det(E) = %e, diag(B). = %e\n', detA, detE, diagB);
        
        tic; % time measure
        B = -B;
        E = -E;
        p = size(A,1);
        q = size(B,1);
        R = zeros(p,q);
        L = zeros(p,q);
        % STEP 1: Transform (A,D) and (B,E) into simpler form. And this has
        % been already done for the input.
        
        for i = 1:q
        % STEP 2: Solve for ri using backward substitution.
        rhs1 = zeros(p,1);
        rhs2 = zeros(p,1);
        if( i > 1 )
            rhs1 = sum(B(1:i-1,i).'.*L(:,1:i-1),2);
            rhs2 = sum(E(1:i-1,i).'.*L(:,1:i-1),2);
        end
        b = C(:,i) - rhs1;
        R(:,i) = backSub(A,b);
        % STEP 3: Compute li and repeat STEP 2 and 3 from i = 1 to q
        L(:,i) = -1/E(i,i)*(-F(:,i)+D*R(:,i)+rhs2);
        end
        time = toc;
        % STEP 3-(2): Check your soln.
        checkC = A*R + L*B;
        checkF = D*R + L*E;
        
%         cerror = sum(sum(abs((checkC - C)./C)));
%         ferror = sum(sum(abs((checkF - F)./F)));
        cerror = sum(sum(abs((checkC - C))));
        ferror = sum(sum(abs((checkF - F))));
        fprintf('Time for solving (GS-algorithm): %f [sec]\n', time);
        tol = 1e-6;
        if((abs(cerror) > tol)||(abs(ferror) > tol))
            fprintf('GS alogrithm\n');
            fprintf('C error: %f\nF error: %f\n', cerror, ferror);
        end
end % End of the function