function [x,u,du,d2u] = genBlasius
xa = 0; xb = 7;
N = 5000;

xmesh = linspace(xa,xb,N);
solinit = bvpinit(xmesh, @guess);

sol = bvp4c(@bvpfcn, @bcfcn, solinit);

x = solinit.x;
u = sol.y(2,:);
du = sol.y(3,:);
d2u = sol.yp(3,:);

end

function dydx = bvpfcn(x,y)
% dydx = zeros(1,3);
dydx = [y(2)
        y(3)
       -y(1)*y(3)/2];
end

function res = bcfcn(ya,yb)
res = [ya(1)
       ya(2)
       yb(2)-1];
end

function g = guess(x)
g = [0
     0
     0];
end