function a = intFcn(x,n, time, axial, azim, N1, rad)
% time, axial, azim, N1, rad
% n = 0: volume integral
% n = 1: surface integral
% Make sure to match cutoff and cutoff factor globally.
cutoff = 0.001;
cutofffactor = 1;
if ( n == 0 )
    [~, x] = ifft3DIntPart_original(x, time, axial, azim, N1, 0, 1, rad, 1, 1, cutoff, cutofffactor);
    [~, a] = ifft3DIntPart_original(x, time, axial, azim, N1, 0, 1, rad, 3, 0, cutoff, cutofffactor);
elseif ( n == 1)
    [~, x] = ifft3DIntPart_original(x, time, axial, azim, N1, 0, 1, rad, 1, 0, cutoff, cutofffactor);
    [~, a] = ifft3DIntPart_original(x, time, axial, azim, N1, 0, 1, rad, 3, 0, cutoff, cutofffactor);
end
end