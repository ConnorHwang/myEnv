function comp = compactness(eval,evec,rp1)
global L1 N1
intpres = 200;
intpx = linspace(0,L1,intpres);

evec = evec(:); % Make sure evec is a columm vector
% Integration of the magnitude over the domain
chi = evec .* conj(evec);
intx = interp1(rp1,chi(1:N1),intpx,'pchip');
inty = interp1(rp1,chi(N1+1:2*N1),intpx,'pchip');
intz = interp1(rp1,chi(2*N1+1:3*N1),intpx,'pchip');
chihat  = trapz(intpx,intx) + trapz(intpx,inty) + trapz(intpx,intz);
% chihat  = trapz(intpx,inty);
evec_norm = evec / sqrt(chihat) * sqrt(eval);
% evec_norm = evec / sqrt(chihat);

% Integration of the mode, gamma
% gamx = interp1(rp1,evec(1:N1),intpx,'pchip');
gamy = interp1(rp1,evec(N1+1:2*N1),intpx,'pchip');
% gamz = interp1(rp1,evec(2*N1+1:3*N1),intpx,'pchip');
% gamma  = trapz(intpx,gamx) + trapz(intpx,gamy) + trapz(intpx,gamz);
gamma  = trapz(intpx,gamy);
% scale = conj(gamma) / abs(gamma);
scale = gamma / abs(gamma);
% scale = 1;

% output of the evec
comp = evec_norm * scale;

% Normalization, maximum value of 'u' to be 1.
% normal = max(real(comp(1:N1)));
% comp = comp / normal / 16;
end