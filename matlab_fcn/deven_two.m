function d1 = deven_two(N,L)

num = round(abs(N));
d1 = zeros(num,num);

for i = 0:(num-1)
    for j = (i+1):2:(num-1)
        d1(i+1,j+1) = 2*j*(2/L); % We don't need Linear Scaling HERE
    end;
end;

d1(1,:) = d1(1,:)/2;

end