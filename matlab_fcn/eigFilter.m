function [eval_,evec_] = eigFilter(eval,evec,tol)
% Find indices
nanind = isnan(eval);
fprintf('NaN: %d\n', sum(nanind));
infind = isinf(eval);
fprintf('inf: %d\n', sum(infind));
% zero-out
logical = nanind + infind;
eval = eval(~logical);
evec = evec(:,~logical);
% Any filter restriction
if(tol~=0)
    logical2 = abs(eval) > tol;
    eval = eval(~logical2);
    evec = evec(:,~logical2);
end
% output
eval_ = eval;
evec_ = evec;
end