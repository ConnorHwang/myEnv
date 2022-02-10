
%  The script file orrsom.m computes the eigenvalues of the Orr-Sommerfeld
%  equation using NxN Chebyshev differentiation matrices.

% S.C. Reddy, J.A.C. Weideman 1998.  Code modified to display output
% by JACW, May 2003.

% Modified by Hanul Connor Hwang on Feb. 4. 2020.
clearvars; close all; clc
% N = input(' Order of the differentiation matrix: N = ? ');
% R = input(' Reynolds number: R = ? ');
N = 60;
R = 2000;
i = sqrt(-1);

[x,DM] = chebdif(N+2,2);                       % Compute second derivative
    D2 = DM(2:N+1,2:N+1,2);                    % Enforce Dirichlet BCs
                                     
[x,D4] = cheb4c(N+2);                          % Compute fourth derivative
     I = eye(size(D4));                        % Identity matrix

A = (D4-2*D2+I)/R-2*i*I-i*diag(1-x.^2)*(D2-I); % Set up A and B matrices
B = D2-I;

e1 = eig(A,B);                                  % Compute eigenvalues

% e1 = diag(e1);
[~,id] = sort(-real(e1)); % assuming real eigenvalues
e1 = e1(id);
% [m,l] = max(real(e1));                          % Find eigenvalue of largest
% disp('Eigenvalue with largest real part = ')   % real part
% disp(e1(l))

%% What if we explicitly add boundary condition?
% By adding zero rows and explicity B.Cs
N = 60;
R = 2000;
i = sqrt(-1);

[x,DM] = chebdif(N+2,4);   % Compute second derivativee
D1 = DM(:,:,1);    
D2 = DM(:,:,2);                    
D4 = DM(:,:,4);                    
I = eye(size(D4));                        % Identity matrix
     
A = (D4-2*D2+I)/R-2*i*I-i*diag(1-x.^2)*(D2-I); % Set up A and B matrices
B = D2-I;

% B.Cs
B(1,:) = 0.;
B(1,:) = [1 zeros(1,N+1)];
B(N+2,:) = 0.;
B(N+2,:) = [zeros(1,N+1) 1];
B(2,:) = 0.;
B(2,:) = D1(1,:);
B(N+1,:) = 0.;
B(N+1,:) = D1(N+2,:);
A(1,:) = 0.;
A(N+2,:) = 0.;
A(2,:) = 0.;
A(N+1,:) = 0.;

e2 = eig(A,B);                                  % Compute eigenvalues
[~,id] = sort(-real(e2)); % assuming real eigenvalues
e2 = e2(id);

% [m,l] = max(real(e2));                          % Find eigenvalue of largest
% disp('Eigenvalue with largest real part = ')   % real part
% disp(e2(l))

%% compare
csize = max(size(e1,1),size(e2,1));
comp = zeros(csize,2);
comp(1:size(e1,1),1) = e1;
comp(1:size(e2,1),2) = e2;