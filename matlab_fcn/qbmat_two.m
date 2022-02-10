function [qb,invF] = qbmat_two(M,xu,e)

%Phase 1

work = M*xu;
A = work'*work;

for i = 1:size(A,1)
    
    xu(:,i) = xu(:,i)/sqrt(A(i,i));
    
end

work = M*xu;
A = work'*work;

%Phase2
F = chol(A);
invF = inv(F);

%Phase 3: 

qb = -sqrt(-1)*F*diag(e)/F;

end 