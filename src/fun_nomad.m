function out = fun_nomad(x, A, lb, ub)

c1 = A*x - ub;
c2 = lb - A*x;

out = [fun(x), c1', c2']; 