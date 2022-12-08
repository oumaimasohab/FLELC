% Find the projection of $-g$ on the subspace defined by the positive base
function projg = projected_gradient(g, D, Z, H)
options = optimset('Display', 'off', 'Algorithm', 'active-set','TolFun', 2*eps);
[coef, fval]  = quadprog(D'*D, (H*Z'*g)'*D,[],[],[],[],zeros(size(D,2),1),Inf*ones(size(D,2),1), zeros(size(D,2),1),options);
projg = D*coef;
