function b = isfeasible(x,presbounds, preslineq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas)

if presbounds
    xl = (min(x-lb) > -tol_feas);
    xu = (max(x-ub) < tol_feas);
else
    xl = 1;
    xu = 1;
end
%
if ~isempty(Aeq)
    Ax = (norm(Aeq*x-beq,Inf) < tol_feas);
else
    Ax = 1;
end
%
if preslineq
    xli = (min(Aineq*x-lbineq) > -tol_feas);
    xui = (max(Aineq*x-ubineq) < tol_feas);
    Aix = (xli && xui);
else
    Aix = 1;
end
b =  xl && xu && Ax && Aix;