function [x0, lb, ub, Aeq, beq, Aineq, lbineq, ubineq,  A, d, e, bbtype] = get_problem(prob)

[x0,lb,ub,~,d,e]=spamfunc(['Cute_nl/' prob]);

[~,A]=spamfunc(x0,1);
A = full(A);

Aineq = [];
Aeq = [];
lbineq =[];
ubineq = [];
beq = [];
bbtype = ['OBJ'];
to_remove = [];
for j=1:size(A, 1)
    bbtype = [bbtype ' PB PB'];

    if (d(j) ~= e(j))
        [t, ind, equal] = is_basis(A(j, :));
        if (t)
            if (equal)
                lower = d(j);
                upper = e(j);
            else
                lower = -e(j);
                upper = -d(j);
            end
            lb(ind) = max(lower, lb(ind));
            ub(ind) = min(upper, ub(ind));
            to_remove = [to_remove, j];
        else
            Aineq = [Aineq; A(j,:)];
            lbineq = [lbineq; d(j)];
            ubineq = [ubineq; e(j)];
        end

    else
        Aeq = [Aeq; A(j, :)];
        beq = [beq; d(j)];
    end
end


A(to_remove, :) = [];
d(to_remove) = [];
e(to_remove) = [];
        


x0 = project_point(x0,lb,ub,[], [],A,d,e);

end

function [t, i, equal] = is_basis(a)

posa = abs(a);

ind = find(posa == 1);

if (length(ind) == 1 && sum(posa) == 1)
    t = true;
    i = ind;

else
    t = false;
    i = -1;
end

equal = all(a == posa);

end

