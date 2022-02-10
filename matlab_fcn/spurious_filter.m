 function [evl_new, evc_new] = spurious_filter(evl, evc, tol, alp, u)

% Filters the eigenvalues
% tol is tolerence for the filter
% alp is streamwise wave number
% u is phase1(under domain) velocity vector

N = max(size(evl,1), size(evl,2));
% fprintf("*** Spurious mode filtering ...\n");
% fprintf("Original number of eigenvalues is %d\n", N);
evl_new = [];
evc_new = [];
j = 1;
N_inf = 0;
for i = 1:N
   if( abs(evl(i)) < tol )
       continue;
   % For the spurious mode that coming from the kinematic condition
   % Considering both cases for 'omega' and 'alp'
   % For temporal stability exclude 'alp * u1(1)', for kinematic cond.
   % For spatial stability exclude 'omega / u1(1)', for kinematic cond.
%    elseif( abs(evl(i)) == abs(alp*u(1)) || abs(evl(i)) == abs(u(1)) ...
%         || abs(evl(i)) == abs(alp/u(1))                             )
%        continue;
   % get rid of the interface mode
   elseif( abs(imag(evl(i))) < 10^-10 )
       continue;
   % For too many chebyshev modes, excludes abnormally large eigenvalues
   elseif( abs(evl(i)) > 10^3 )
       N_inf = N_inf + 1;
       continue;
   elseif( abs(evl(i)) == inf )
       N_inf = N_inf + 1;
   else
       evl_new(j) = evl(i);
       evc_new(:,j) = evc(:,i);
       j = j+1;
   end
end

evl_new = evl_new.';

N_after = max(size(evl_new,1), size(evl_new,2));

% fprintf("Number of infinite eigenvalues is %d \n", N_inf);
% fprintf("Number of filtered eigenvalues is %d \n", N - N_after - N_inf);
% fprintf("Number of eigenvalues after filtering is %d\n", N_after);

end