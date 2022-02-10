function [D0m,D1m,D2m,D3m,D4m] = DmatC(N,L)

% initialize

num = (round(abs(N))-1);

% y coordinate in Chebyshev space

vec = (0:1:num)';
y_cheb = cos(pi*vec/num);

% y coordinate in physical space

y_phy = L*(y_cheb + 1)/2;

% Derivative of y

dy_cheb_dy_phy = 2/L;

% Create D0

D0 = [];
for j=0:1:num
    D0 = [D0 cos(j*pi*vec/num)];
end;

% Create higher derivative matrices

lv = length(vec);
D1 = [zeros(lv,1) D0(:,1) 4*D0(:,2)];
D2 = [zeros(lv,1) zeros(lv,1) 4*D0(:,1)];
D3 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
D4 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
for j=3:num

    D1 = [D1 2*j*D0(:,j)+j*D1(:,j-1)/(j-2)];
    D2 = [D2 2*j*D1(:,j)+j*D2(:,j-1)/(j-2)];
    D3 = [D3 2*j*D2(:,j)+j*D3(:,j-1)/(j-2)];
    D4 = [D4 2*j*D3(:,j)+j*D4(:,j-1)/(j-2)];

end;

D0m = D0;
D1m = D1*(dy_cheb_dy_phy);
D2m = D2*(dy_cheb_dy_phy)^2;
D3m = D3*(dy_cheb_dy_phy)^3;
D4m = D4*(dy_cheb_dy_phy)^4;

end