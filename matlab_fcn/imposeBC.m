function B = imposeBC(A,N,brow,bcol,erow,bc,drow)
% Do not use this function on the block matrix
global D01 D11
er = 200;
% impose BC
for i = 1:size(brow,2)
    for j = 1:size(bcol,2)
        col = (bcol(j)-1)*N+1:bcol(j)*N;
        for k = 1:size(erow,2)
            row = N*(brow(i)-1)+mod(erow(k),N+1);
            switch(bc)
                case 0
                    if(drow > 0)
                        b = D01(1,:);
                    else
                        b = D01(end,:);
                    end % End of if
                case 1
                    if(drow > 0)
                        b = D11(1,:);
                    else
                        b = D11(end,:);
                    end % End of if
                case 2
                    if(drow > 0)
                        b = er*D01(1,:);
                    else
                        b = er*D01(end,:);
                    end % End of if
                case 3
                    if(drow > 0)
                        b = er*D11(1,:);
                    else
                        b = er*D11(end,:);
                    end % End of if
                otherwise
                    b = zeros(1,N);
            end % End of 'switch'
            % zero-out the entire row
            A(row,:)   = 0;
            A(row,col) = b;
        end
    end
end
B = A;
end