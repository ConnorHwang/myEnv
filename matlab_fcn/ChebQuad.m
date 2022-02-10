function wi = ChebQuad(N1,L1)
% 'ChebQuad' returns coefficient 'wi' that sum of wi*fi would give the
% integration of 'f' over the domain of size L1 with N1 Chebyshev points.
j = 0:N1-1; j = j(:);
xch = cos(pi*(j)/(N1-1)); xch = xch(:);
wi = (pi/(N1-1)*(sin(j/(N1-1)*pi).^2))./(sqrt(1-xch.^2));
wi = wi(:)*L1/2;
er = 1e-8;
wi(isinf(wi))=er; wi(isnan(wi))=er;
end