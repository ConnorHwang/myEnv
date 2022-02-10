function [f,df,h] = genFilt(rp1,a,del99,option)
% set parameters
switch(option)
    case 1
        y1 = 0; y2 = 5;
    case 2
        y1 = 5; y2 = 10;
    case 3
        y1 = 10; y2 = 15;
    case 4
        y1 = 15; y2 = 20;
end % End of 'switch'
% function h
h = ones(size(rp1,1),1);
f = 1/pi*(atan(a*(rp1-y1))-atan(a*(rp1-y2)));
df = a/pi*(1./((a*(rp1-y1)).^2+1)-1./((a*(rp1-y2)).^2+1));
end