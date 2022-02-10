function B = rmRowCol(A,rmind)
% remove rows and columns
% Input -------------------------------------------------------------------
% A    : Original matrix
% rmind: Indices that you want to remove
% Output ------------------------------------------------------------------
% B    : Removed matrix
%--------------------------------------------------------------------------
% remove rows
A(rmind,:) = [];
% remove cols
A(:,rmind) = [];
B = A; % return A
end % End of the function