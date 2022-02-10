function [U,dxU,dyU,dyyU,V,dyV,del99] = genBaseFlowGlobal(rgrid,xgrid,sscale)
% When 'eta' used as a characteristic length, predicted normal velocity for an infinite value.
vinf = 0.8604*sscale;
Nx   = max(size(xgrid,1),size(xgrid,2));
Ny   = max(size(rgrid,1),size(rgrid,2));
U    = zeros(Ny,Nx);
dxU  = zeros(Ny,Nx);
dyU  = zeros(Ny,Nx);
dyyU = zeros(Ny,Nx);
V    = zeros(Ny,Nx);
dyV  = zeros(Ny,Nx);
%%{
% Determine how the Blasius B.L. grows
del_ini = 3.5*sqrt(2)*sscale; % 3.5 is del99 when sqrt(2)*delta_0 as a characteristics length
del_fin = 7.8*sqrt(2)*sscale;
xmax = xgrid(1);
x1 = xmax/((del_fin/del_ini)^2-1);
del99 = del_ini*sqrt(1+xgrid/x1);

% The data uses delta_0, so f''(0) ~ 0.332
load BlasiusProfile.mat zBlas fBlas

vBlas = flipud((zBlas.*fBlas(:,2)-fBlas(:,1))/2);
% vBlas = flipud((zBlas.*fBlas(:,2)-fBlas(:,1)));
zBlas = flipud(zBlas);
fBlas = flipud(fBlas);
dvBlas = gradient(vBlas,zBlas);
% figure;plot(vBlas,zBlas,dvBlas,zBlas);
% at x=0
for i = 1:Nx
%--------------------------------------------------------------------------
% Commented out for parallel flow test
if(Nx==1)
    C = 5/3.5/sqrt(2); 
%     C = 1;
%--------------------------------------------------------------------------
% Comment for parallel flow test
else
    C = del99(i)/3.5/sqrt(2);
    v = interp1(zBlas,vBlas,rgrid/C,'linear',vinf); % Uses delta_0, 
    dv = (1/C)*interp1(zBlas,dvBlas,rgrid/C,'linear',0); 
    V(:,i) = v;
    dyV(:,i) = dv; 
end
%--------------------------------------------------------------------------
    u = interp1(zBlas,fBlas(:,2),rgrid/C,'linear',1);
    du = (1/C)*interp1(zBlas,fBlas(:,3),rgrid/C,'linear',0); % Note that f''(eta = 0) = 0.3321
    d2u = (1/C)^2*interp1(zBlas,fBlas(:,4),rgrid/C,'linear',0); d2u(1) = 0;
    d2u(1) = 0;
    U(:,i) = u;
    dyU(:,i) = du;
    dyyU(:,i) = d2u;
%--------------------------------------------------------------------------
end
%}
% dxU
%{
Ninterp = 1000;
Uinterp = zeros(Ny,Ninterp);
dxUinterp = zeros(Ny,Ninterp);
xinterp = linspace(xgrid(1),xgrid(end),Ninterp);
del99_interp = del_ini*sqrt(1+xinterp/x1);
for i = 1:Ninterp
%     C = del99(i)/3.5;
    C = del99_interp(i)/3.5/sqrt(2);
    u = interp1(zBlas,fBlas(:,2),rgrid/C,'linear',1);
    Uinterp(:,i) = u;
end
% interpolation
for i = 1:Ny
    dxUinterp(i,:) = gradient(Uinterp(i,:).',xinterp.');
    dxU(i,:) = interp1(xinterp.',dxUinterp(i,:).',xgrid,'linear').';
end
%}
%%
% ch = 9;
% figure; plot(xinterp,Uinterp(ch,:),'k-'); hold on;
% plot(xinterp,dxUinterp(ch,:),'b-');
% a = Uinterp(ch,:);
% b = gradient(a,xinterp);
% plot(xinterp,b)
%% plot baseflow
%{
figure;
pcolor(xgrid,rgrid,U); hold on;
title('U(x,y)');
shading interp
colorbar
colormap jet
plot(xgrid,del99,'k-','linew',2);
axis square
box on;
% Base flow, V(x,y)
figure;
pcolor(xgrid,rgrid,dxU); hold on;
title('V(x,y)');
shading interp
colorbar
colormap jet
plot(xgrid,del99,'k-','linew',2);
axis square
box on;
%}
end % End of the function