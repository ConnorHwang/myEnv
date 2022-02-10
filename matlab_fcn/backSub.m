function x = backSub(U,b)
S=size(U);
m=S(1);
if S(1)~=S(2)
    error('matrix mast be square')
end
% Make sure input 'b' is a vector form
b = b(:);

% First value
x = zeros(m,1);
x(m,1)=b(m)/U(m,m);

% Bacward substitution
for k=m-1:-1:1  
        x1=(b(k)-sum(U(k,k+1:m)*x(k+1:m)))/U(k,k);
        x(k)=x1;
end
end