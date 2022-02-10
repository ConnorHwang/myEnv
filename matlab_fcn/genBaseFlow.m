function [u,du,d2u] = genBaseFlow(r,sscale,option)
switch(option)
    case 1 % Blasius BL
        % Note that eta(2*del_0) = 3.5, f'(eta) ~ 0.99, eta = 5, f'(eta) = 0.99994;
        [bu,u,du] = BlasiusSol(r,sscale,'no');
        d2u = -du.*bu;
    case 2
        % The data uses delta_0, so f''(0) ~ 0.332
        load BlasiusProfile.mat zBlas fBlas
        % Interpolate
        zBlas = flipud(zBlas);
        fBlas = flipud(fBlas);
        u = interp1(zBlas,fBlas(:,2),sscale*r,'linear',1);
        du = interp1(zBlas,fBlas(:,3),sscale*r,'linear',0); % Note that f''(eta = 0) = 0.3321
        d2u = interp1(zBlas,fBlas(:,4),sscale*r,'linear',0);
        d2u(1) = 0;
        du = sscale*du;
        d2u = sscale^2*d2u;
%         du = sqrt(2)*du;
%         d2u = 2*d2u;
end % End of 'switch'
% % Plotting profile from the loaded data
% figure;
% %         plot(fBlas(:,2),zBlas,fBlas(:,3),zBlas,fBlas(:,4),zBlas);
% plot(u,r,'*-',du,r,'*-',d2u,r,'*-');
% legend('u','du','d2u');
% box on; grid on;
% set(gca,'linew',2,'FontSize',15);
end