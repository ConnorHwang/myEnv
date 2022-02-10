function saveShiftedInit(N2,L1,L2,r1,r2,u1,u2,Uini,Vini,Wini,Fini,Nr,output,path,do_save)
% Save initial optimized initial perturbations that are phase
% shifted for 'f' to have zero magnitude in physical space at t = 0.
% Output: 
% Define interpolation values
r_tot = [r2(1:end-1); r1];
u_tot = [u2(1:end-1); u1];
rEdge = L1 + L2;
r     = (0:Nr-1)*rEdge/Nr; r = r(:);
input_base = interp1(r_tot, u_tot, r, 'pchip');
% For loop for each theta value
Urbuf = real(Uini);
Uibuf = imag(Uini);
Vrbuf = real(Vini);
Vibuf = imag(Vini);
Wrbuf = real(Wini);
Wibuf = imag(Wini);
% Eliminate the overlapping region for interpolation
Urbuf(N2) = [];
Uibuf(N2) = [];
Vrbuf(N2) = [];
Vibuf(N2) = [];
Wrbuf(N2) = [];
Wibuf(N2) = [];
% Interpolate
MUr = interp1(r_tot, Urbuf, r, 'pchip');
MUi = interp1(r_tot, Uibuf, r, 'pchip');
MVr = interp1(r_tot, Vrbuf, r, 'pchip');
MVi = interp1(r_tot, Vibuf, r, 'pchip');
MWr = interp1(r_tot, Wrbuf, r, 'pchip');
MWi = interp1(r_tot, Wibuf, r, 'pchip');
MB = input_base;
    
fprintf("f(t = 0) = %f + %f(i) \n", real(Fini), imag(Fini));
fprintf("|f|      = %f \n", abs(Fini));

% Plot figure to compare with SVD results.
figure; hold on;
plot(input_base,r);
plot(MUr,r);
plot(MUi,r);
plot(MVr,r);
plot(MVi,r);
plot(MWr,r);
plot(MWi,r);
legend('base','ur','ui','vr','vi','wr','wi');
% figure;
% plot(abs(blkdiag(D0_2,D0_1)*iniu([U2vec U1vec])),[r2; r1],'k','linew',2); hold on; grid on; box on;
% plot(abs(blkdiag(D0_2,D0_1)*iniu([V2vec V1vec])),[r2; r1],'k--','linew',2);
% plot(abs(blkdiag(D0_2,D0_1)*iniu([W2vec W1vec])),[r2; r1],'k-.','linew',2);

if(do_save == 1)
    % Save file
    dlmwrite(fullfile(path,[output '_base.txt']),MB, 'delimiter','\n');
    dlmwrite(fullfile(path,[output '_ur.txt']),  MUr,'delimiter','\n','newline','pc');
    dlmwrite(fullfile(path,[output '_ui.txt']),  MUi,'delimiter','\n','newline','pc');
    dlmwrite(fullfile(path,[output '_vr.txt']),  MVr,'delimiter','\n','newline','pc');
    dlmwrite(fullfile(path,[output '_vi.txt']),  MVi,'delimiter','\n','newline','pc');
    dlmwrite(fullfile(path,[output '_wr.txt']),  MWr,'delimiter','\n','newline','pc');
    dlmwrite(fullfile(path,[output '_wi.txt']),  MWi,'delimiter','\n','newline','pc');
    fprintf("Simulation input file, %s, saved!\n", output);
elseif(do_save == 2) % Save for Cascade code
    R = 0.0039;
    MB((abs(MB) < 1e-4)) = 0;
    dlmwrite(fullfile(path,'umean.dat'),[r*R,MB] ,'delimiter','\t');
    dlmwrite(fullfile(path,'ur.dat'),  [r*R,MUr],'delimiter','\t');
    dlmwrite(fullfile(path,'ui.dat'),  [r*R,MUi],'delimiter','\t');
    dlmwrite(fullfile(path,'vr.dat'),  [r*R,MVr],'delimiter','\t');
    dlmwrite(fullfile(path,'vi.dat'),  [r*R,MVi],'delimiter','\t');
    dlmwrite(fullfile(path,'wr.dat'),  [r*R,MWr],'delimiter','\t');
    dlmwrite(fullfile(path,'wi.dat'),  [r*R,MWi],'delimiter','\t');
    fprintf("Simulation input file, %s, saved!\n", output);
end
end




