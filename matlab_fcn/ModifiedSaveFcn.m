function ModifiedSaveFcn(N2,L1,L2,r1,r2,n,u1,u2,Uini,Vini,Wini,Fini,Nr,Nthe,output)
% Define interpolation values
kappa   = 1;% baseflow scalar
epsilon = 1e-3; % perturbation amp. scalar
r_tot = [r2(1:end-1); r1];
u_tot = [u2(1:end-1); u1];
theta = (0:Nthe-1)*2*pi/Nthe;
rEdge = L1 + L2;
r     = (0:Nr-1)*rEdge/Nr;
input_base = interp1(r_tot, u_tot, r, 'pchip');
% Container for initial velocity
MB = zeros(Nthe*Nr,1); % baseflow without disturbances
MU = zeros(Nthe*Nr,1);
MV = zeros(Nthe*Nr,1);
MW = zeros(Nthe*Nr,1);
MF = zeros(Nthe   ,1);
% For loop for each theta value
for i = 1:Nthe
    Ubuf = Uini*sin(n*theta(i));
    Vbuf = Vini*sin(n*theta(i));
    Wbuf = Wini*sin(n*theta(i));
    Fbuf = Fini*sin(n*theta(i));
    Ubuf(N2) = [];
    Vbuf(N2) = [];
    Wbuf(N2) = [];
    % Interpolate
    input_u = interp1(r_tot, Ubuf, r, 'pchip');
    input_v = interp1(r_tot, Vbuf, r, 'pchip');
    input_w = interp1(r_tot, Wbuf, r, 'pchip');
    input_f = Fbuf;

    MB((i-1)*Nr+1:(i-1)*Nr+Nr,:) = kappa*input_base;
    MU((i-1)*Nr+1:(i-1)*Nr+Nr,:) = input_u*epsilon;
    MV((i-1)*Nr+1:(i-1)*Nr+Nr,:) = input_v*epsilon;
    MW((i-1)*Nr+1:(i-1)*Nr+Nr,:) = input_w*epsilon;
    MF(i                       ) = input_f*epsilon;
end
fprintf("method(b), Famp = %e\n", max(MF));
% Plot figure to compare with SVD results.
% figure;
% plot(abs(blkdiag(D0_2,D0_1)*iniu([U2vec U1vec])),[r2; r1],'k','linew',2); hold on; grid on; box on;
% plot(abs(blkdiag(D0_2,D0_1)*iniu([V2vec V1vec])),[r2; r1],'k--','linew',2);
% plot(abs(blkdiag(D0_2,D0_1)*iniu([W2vec W1vec])),[r2; r1],'k-.','linew',2);
% Save file
dlmwrite([output '_base.txt'],MB,'newline','pc');
dlmwrite([output '_u.txt'],MU,'newline','pc');
dlmwrite([output '_v.txt'],MV,'newline','pc');
dlmwrite([output '_w.txt'],MW,'newline','pc');
dlmwrite([output '_f.txt'],MF,'newline','pc');
fprintf("Simulation input file, %s, saved!\n", output);
end




