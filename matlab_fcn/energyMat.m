function M = energyMat(N,L,option)
% case1: Schmid's D matrices, case2: Weideman's.
% State-space energy matrix.
% input: state-space vector, [v_hat eta_hat]^T.
% output: Domain integration of the kinetic energy. \int_domain
% 1/2*(u_hat_i^2) dy.
global kx kz D11
ak2 = kx^2+kz^2;
% M = eye(N+N,N+N);
switch(option)
    case 1
        Cos = intMat(N,0,L,'cart');
        Dos = derCheb(N,L);
        Wos = (Dos'*Cos*Dos+ak2*Cos)/ak2;
        Wsq = intMat(N,0,L,'cart')/ak2;
        % Wos = (Dos'*Cos*Dos+ak2*Cos)/(ak2);
        % Wsq = two_two(N)/ak2;
        [u,s,~]=svd(Wos); s=sqrt(diag(s));
        Mos=diag(s)*u';
        [u,s,~]=svd(Wsq); s=sqrt(diag(s));
        Msq = diag(s)*u';
        M = [Mos zeros(N,N); zeros(N,N) Msq];
    case 2
        Cos = ChebQuad(N,L); Cos = Cos(2:N-1);
        Dos = D11;
        Wos = (Dos'*diag(Cos)*Dos+ak2*diag(Cos))/ak2;
        Wsq = diag(Cos)/ak2;
        % Wos = (Dos'*Cos*Dos+ak2*Cos)/(ak2);
        % Wsq = two_two(N)/ak2;
        [u,s,~]=svd(Wos); s=sqrt(diag(s));
        Mos=diag(s)*u';
        [u,s,~]=svd(Wsq); s=sqrt(diag(s));
        Msq = diag(s)*u';
        M = [Mos zeros(N-2); zeros(N-2) Msq];
    case 3
        Cos = ChebQuad(N,L);
        Dos = D11;
        Wos = (Dos'*diag(Cos)*Dos+ak2*diag(Cos))/ak2;
        Wsq = diag(Cos)/ak2;
        [u,s,~]=svd(Wos); s=sqrt(diag(s));
        Mos=diag(s)*u';
        [u,s,~]=svd(Wsq); s=sqrt(diag(s));
        Msq = diag(s)*u';
        M = [Mos zeros(N); zeros(N) Msq];
end
end