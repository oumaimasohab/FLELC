function [x_proj] = project_point(x, lb,ub,Aeq,beq,Aineq,lbineq,ubineq)


A = [Aineq; -Aineq];

b = [ubineq; -lbineq];

x_proj = quadprog(eye(length(x)),-x,A,b,Aeq,beq,lb,ub);

end
