function [alpha,fval,func, xi, p] = armijo(fc,f,x,H,g_current,xi,ds_step,presbounds, preslineq, presleq, lb,...
    ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas, Ai,ni, li, ui, Z, W, fle)
%% Compute the backtracking Armijo search for finding alpha
alpha = 1;
k     = 0;
fval = 0;
eta = 1e-8;
c_2 = 1e-10;
tau = 0.5;
func = fc;
% xi = 1;
xi_init = xi;
% fprintf('             Sufficient decrease with backtracking            \n')
% fprintf('iter       alpha      norm-next-step\n') 
% [g, p] = project_on_approx_tangent_cone(x, g_current, H, ds_step, preslineq, presbounds, presleq, Aeq, Ai, ni, li, ui, Z, W);

while (1)
    xi = alpha*xi_init;
%     Project the gradient along with the step p_current
    [g, p] = project_on_approx_tangent_cone(x, g_current, H, xi, preslineq, presbounds, presleq, Aeq, Ai, ni, li, ui, Z, W);
    xn = x + xi * p/norm(p);
    %use of min
    if (fle == 1 && xi < min(1e-2, 1e-3*ds_step))
%         disp("Armijo failed:")
%         xi
%         disp("alpha:")
%         alpha
        alpha = 0;
        func  = fc;
        xi = xi_init/2;
        return
    end
    if (xi <= 1e-16)
        alpha = 0;
        xi = xi_init/2;
        func  = fc;
        return
    end

    if  isfeasible(xn,presbounds, preslineq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas)
        try
            func = feval(f, xn);
        catch
            alpha = 0;
            warning(' SD has trouble computing f. exiting program.');
            disp(norm(p))
            disp(norm(g))
            fprintf('iter       alpha   norm step   norm-next-step\n') 
            fprintf('%5g %14.7e %14.7e %14.7e  \n',k, alpha,norm(x), norm(x+alpha*p))
%             pause
            return
        end
    
    %    fprintf('%5g %14.7e %14.7e \n',k, alpha, norm(x+alpha*p))
        fval = fval + 1;

        if (func <= fc - eta*xi*g'*p/norm(p) - c_2* xi^2)
            alpha = xi;
%             Bound constraints
%             xi = min(100,10*xi);
            % Linear equality
            xi = min(10,xi_init*2);
%             xi = xi_init;
        return
        end
    end

    alpha = tau * alpha;
    k     = k +1;
end