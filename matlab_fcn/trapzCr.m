function cr = trapzCr(phi,L1,N1,rp1,del99)
[evec,eval] = eig(phi);
eval = real(diag(eval));
tt = 0;
intpres = 1000;
for i = 1:size(evec,2)
vel_til = evec(:,i);
u_til = vel_til(1:N1);
v_til = vel_til(N1+1:2*N1);
w_til = vel_til(2*N1+1:3*N1);
% Integrate up to the boundary layer
intpx = linspace(0,L1,intpres);
uhat = interp1(rp1,abs(sqrt(eval(i))*u_til),intpx,'pchip');
vhat = interp1(rp1,abs(sqrt(eval(i))*v_til),intpx,'pchip');
what = interp1(rp1,abs(sqrt(eval(i))*w_til),intpx,'pchip');
% figure; plot(intpx,uhat,intpx,vhat,intpx,what);
intpart = intpx(intpx<del99);
usum = trapz(intpart,uhat(intpx<del99).^2);
vsum = trapz(intpart,vhat(intpx<del99).^2);
wsum = trapz(intpart,what(intpx<del99).^2);
tt = tt+(usum+vsum+wsum)/2;
end
% Dg_ = diag(rp1<del99);
% Dg = blkdiag(Dg_,Dg_,Dg_);
% nomi = trace(Dg*phi*Dg')/2;
% use chebint?
% fprintf('%f =? %f \n', tt, nomi);
cr = log10(tt);
end % End of the function