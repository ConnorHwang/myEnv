function x=nlize_two(N,x,M)

nc = size(x); nc = nc(2);
tp = x(2*N+1,:);
x(2*N+1,:) = [];

for i = 1:nc
    x(:,i) = x(:,i) / norm(M*x(:,i));
end;

x = [x(1:2*N,:); tp; x(2*N+1:4*N,:)];

end