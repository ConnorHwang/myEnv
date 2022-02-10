function c = two_two(N,L)

num = round(abs(N));
c = zeros(num,num);

for i = (0:num-1)
    for j = (0:num-1)
        if rem(i+j,2) == 0
            p = (1/(1-(i+j)^2)+1/(1-(i-j)^2))*(L/2); % 1/2 due to the linear mapping
            c(i+1,j+1)=p;
        else
            c(i+1,j+1)=0;
        end;
    end;
end;

end