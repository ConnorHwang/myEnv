function [n1, n2] = n9_two(e,a)

n1 = 1;

while imag(e(n1))>.5;
    n1 = n1+1;
end;

n2 = n1;

while imag(e(n2))>a,
    n2 = n2+1;
    if n2 == size(e,2)
        n2 = n2+1;
        break;
    end
end;

n2 = n2-1;

end