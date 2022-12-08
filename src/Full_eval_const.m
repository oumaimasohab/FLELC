function [x_current,f_current,flag, output, typeofeval, list_xi, list_alpha, list_beta] = Full_eval_const(fun,x_initial,lb,ub,...
       Aeq,beq,Aineq,lbineq,ubineq,maxevals,tol_feas)
%
% Full eval algorithm for solving problems with linear constraints


%Initialization
seed = 1;
rand('state',seed);
% global Ai
% global li
% global ui
% global Z
% global W
n          = length(x_initial);
func_eval  = 0;
numIter    = 0;
h          = sqrt(eps);
xi         = 1;
x_current  = x_initial;
f_current  = feval(fun,x_current);

grad_comp = 1;

H          = eye(n);
func_eval  =  1;
fvalues    = [f_current];
nfeval     = [func_eval];
tf         = [fun(x_current, 1)];nfeval, 
nLE        = 0;
gamma      = 2;
theta      = 0.5;
cont       = 1;
In         = eye(n);
alpha      = 1e-10;
tol_alpha  = 1e-10;
bfgsUp     = 0;
numFeaspts = 0;
p0 = log(theta)/log(theta/gamma);
beta      = 1;
typeofeval = [];
fval = 5;
% Minimum number of directions required in the unconstrained case
nbunc = ceil(log(1-log(theta)/log(gamma))/log(2));
%

%%%%%%%%%%%%%%%%%%%%%%%% Determination of the presence of bounds %%%%%%%%%%%%%%
%
lowbnd = find(lb > -1e20);
upbnd = find(ub < 1e20);
presbounds = ~isempty(upbnd) || ~isempty(lowbnd);


%
%%%%%%%%%%%%%%%%%%%%%%%% Checking for equality constraints %%%%%%%%%%%%%%%%%%%%
%
if ~isempty(Aeq)
	presleq = 1;
	for i=1:size(Aeq,1)
		auxa = norm(Aeq(i,:));
		Aeq(i,:) = Aeq(i,:)/auxa;
		beq(i) = beq(i)/auxa;
	end
	Z = null(Aeq);
	m = size(Z,2);
%
    if m==0
%         We are already at the solution, or the problem is infeasible
          cont = 0;
    end
else
	presleq = 0;
	Z = In;
    m = n;
end


%
%%%%%%%%%%%%%%%%%%%% Gathering linear inequalities and bounds %%%%%%%%%%%%%%%%%
%
if ~isempty(Aineq)
	preslineq = 1;
	for i=1:size(Aineq,1)
		auxa = norm(Aineq(i,:));
		Aineq(i,:) = Aineq(i,:)/auxa;
		lbineq(i) = lbineq(i)/auxa;
		ubineq(i) = ubineq(i)/auxa;
	end
	Ai = [Aineq;In];
	li = [lbineq;lb];
	ui = [ubineq;ub];
	ni = length(ui);
else
	preslineq = 0;
	Ai = In;
	li = lb;
	ui = ub;
	ni = n;
end
%
% The projected normal vectors
W = Z*Z'*Ai';
FE_step = 0;
in = 0;
IsFstep = 0;
xi = 1;
Im = eye(m);
H = eye(m);
y_current  = zeros(m,1);

list_xi = [];
list_alpha = [];

list_beta = [];

%
notchanged = 1;

while (cont)



    gp = ones(n,1);
    numIter = numIter + 1;

    
    
%     if (numIter ==1)
%         xi =1;
%     else
%         xi = norm(g_pas);
%     end
    %	Computation of the alpha-active and inactive constraints
%     if (activatemin)
%      xi = max(xi, alpha);
%     end
%
%     xi = min(1, norm(g_current));
%     if (IsFstep)
     [Iap,Iam,nbgen] = activeconsAiZ(x_current,xi,li,ui,Ai,Z);
%     else
%          [Iap,Iam,nbgen] = activeconsAiZ(x_current,min(xi,alpha),li,ui,Ai,Z);
%     end
%
    Iip = setdiff(1:ni,Iap);
    Iim = setdiff(1:ni,Iam);
    Ya  = [];
   [g_current,fe] = grad_finitediff(fun,x_current,f_current,h);
        func_eval  = func_eval + fe;
        FE_step = FE_step + 1;
    
    
         if (FE_step == 1)
              g_pas = g_current;
              y_pas = y_current;
         end
    
         in = in + 1;
        
        if  (FE_step > 1)
%             s = x_current - x_pas;
%             y = g_current - g_pas;
%             s = - beta * H*Z'*g_pas;
            s = y_current - y_pas;
            y = Z'*g_current - Z'*g_pas;
            c = y'*s;
            r = 1/c;
            % 
            if (bfgsUp == 0) && (c > 0)
                 H = (c/(y'*y)) * Im;    
            end 
    %                     c
    %                     pause
            if (c > 1e-10 * norm(y) * norm(s))
                    %Update BFGS
                    H = (Im-r*s*y') * H * (Im-r*y*s') + r*(s*s');      % BFGS  (6.17)
            end
            bfgsUp = bfgsUp + 1; 
        end 
    

        FE_step = FE_step + 1;
        numIter = numIter + 1;

         [beta,fval,f_temp, xi, d] =  armijo(f_current,fun,x_current,H,g_current,xi,alpha, ...
            presbounds, preslineq, presleq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas, Ai, ni, li, ui, Z, W, 1);
        
        func_eval = func_eval + fval;
        if (beta == 0 || norm(d) == 0)
            success = 0;

        else
            xtemp     = x_current + beta * d/norm(d);
            g_pas     = g_current;
            x_pas     = x_current;
            y_pas     = y_current;
            success  = 1;
            x_current = xtemp;
            f_current = f_temp;
            y_current =  y_current + Z'*beta*d/norm(d);
        end
        
         fvalues   = [fvalues, f_current];
         nfeval    = [nfeval, func_eval];
         tf        = [tf, fun(x_current, 1)];
%     end
    cont = ( (func_eval < maxevals) && (numIter < 2000*n) );
end


output.numIter = numIter;
output.funcCount = func_eval;
output.fvalues = fvalues;
output.funcCounts = nfeval;
output.truefunc  = tf;
% 
%
if (alpha >=tol_alpha)
	if (func_eval==1)
%		No initial feasible point was found
		flag = 3;
	else 
%		Budget of function evaluations exceeded
		flag = 2;
	end
else
	flag = 1;
end
output.flag = flag;
