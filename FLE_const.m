function [x_current,f_current,flag, output, typeofeval, list_xi, list_alpha, list_beta, list_succ, list_type] = FLE_const(fun,x_initial,lb,ub,...
       Aeq,beq,Aineq,lbineq,ubineq,maxevals,tol_feas, bfgs_use, activatemin)
%
% Full-Low eval algorithm for solving problems with linear constraints


%Initialization
seed = 1;
rand('state',seed);


n          = length(x_initial);
func_eval  = 0;
numIter    = 0;
fvalues    = [];
nfeval     = [];
h          = sqrt(eps);
xi         = 1;
evl        = "FE";
x_current  = x_initial;
f_current  = feval(fun,x_current);
% [g_current,fe] = grad_finitediff(x_current,f_current,h);

% g_pas      = g_current;
% x_pas      = x_current;
H          = eye(n);
func_eval  =  1;
fvalues    = [fvalues, f_current];
nfeval     = [nfeval, func_eval];
tf         = [fun(x_current, 1)];
nLE        = 0;
gamma      = 2;
theta      = 0.5;
cont       = 1;
In         = eye(n);
alpha      = 1;
tol_alpha  = 1e-10;
bfgsUp     = 0;
numFeaspts = 0;
p0 = log(theta)/log(theta/gamma);
beta      = 1;
typeofeval = [];
fval = 5;
% Minimum number of directions required in the unconstrained case
nbunc = ceil(log(1-log(theta)/log(gamma))/log(2));
list_succ = [];
num_success_iter_le = 0;
num_success_iter_fe = 0;
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

list_xi = [];
list_alpha = [];
list_type = [];

list_beta = [];

%
notchanged = 1;

Im = eye(m);
y_current  = zeros(m,1);


H = Im;


% if isfeasible(x_current,presbounds, preslineq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas)
%     disp('you good')
% else
%     pause
% end

% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% fprintf( '  it    iter_t      neval        fvalue       success     ||g||     alpha      xi    beta  \n');

while (cont)

%     disp("current iteration type:")
%     disp(evl)

    notchanged = 1;
    list_xi = [list_xi, xi];
    list_alpha = [list_alpha, alpha];

    gp = ones(n,1);
    
    D = [];
%     if (evl == "FE")
%         typeofeval = [typeofeval, numIter];
%     end
% 
%     if (evl == "FE")
%          [Iap,Iam,nbgen] = activeconsAiZ(x_current,xi,li,ui,Ai,Z);
%     else
         [Iap,Iam,nbgen] = activeconsAiZ(x_current,alpha,li,ui,Ai,Z);
%     end
%
    Iip = setdiff(1:ni,Iap);
    Iim = setdiff(1:ni,Iam);
    Ya  = [];
    if (nbgen == 0)
        % The case where there are no active bounds - situation similar to linear 
% 		equality-constrained (or unconstrained) case
%     
        if (evl == "LE")
            	D = 2*rand(n,max(1,ceil(nbunc/2)))-1;
				if ~preslineq  %&& ~presbounds
					D = Z'*D;
				else
					D = W'*D;
				end
%
				for i=1:size(D,2)
					D(:,i) = D(:,i)/norm(D(:,i));
				end
%
				if ~preslineq % && ~presbounds
					D = Z*D;
				else
					D = W*D;
				end
%
				D = [D -D];
        end
    else
        if (preslineq || (presleq && presbounds))
            [Ys_r,Yc_r,calldd] = gentech(Aeq,Ai,Z,Iap,Iam);
            if isempty(Ys_r)
                Ys = [];
            else
                Ys = Z*Ys_r;
            end
            if isempty(Yc_r)
                Yc = [];
            else
                Yc = Z*Yc_r;
            end	

% 			sumdd = sumdd + calldd;
%
%			Orthonormalization process
%
			if ~isempty(Ys)
				for i=1:size(Ys,2)
					wsi = Ys(:,i);
					Ys(:,i) = wsi/norm(wsi);
				end
%
				if ~isempty(Yc)
%
					Yc = Yc - Ys*Ys'*Yc;
					for i=1:size(Yc,2)
						wci = Yc(:,i);
						Yc(:,i) = wci/norm(wci);
					end
				end
			end

        else
            Is = intersect(Iip,Iim);
			Icp = setdiff(Iip,Is);
			Icm = setdiff(Iim,Is);
			Ys = W(:,Is);
			Yc = [W(:,Icp) -W(:,Icm)];
        end

        dyy = 2*size(Ys,2)+size(Yc,2);    

    %		
        if (evl == "LE" && dyy>0)
    %

            D = [];
    %
    %		    Subspace vectors
            nbdir = ceil(nbgen*p0);
    %
            if ~isempty(Ys)
                Dsb = 2*rand(size(Ys,2),...
                max(1,ceil(nbunc/2)))-1;
                Dsb = Ys*Dsb;
                for i=1:size(Dsb,2)
                    Dsb(:,i) = ...
                    Dsb(:,i)/norm(...
                    Dsb(:,i));
                end
                D = [Dsb -Dsb];
            end
    %		Cone vectors
            if ~isempty(Yc)
                nyc = size(Yc,2);
                nbdir = ceil(nyc*p0);
                Dco = randcolZsearch(Yc,1:nyc,[],nbdir);
                D = [D Dco];
            end
        
        else
			D = [];
        end

	    if ~isempty(Ya)
			    if dirgen>0
				    ca = size(Ya,2);
				    dda = ceil(ca*p0);
				    Ya = Ya(:,dda);
			    end
			    D = [D Ya];
	    end
    end

    if (evl == "FE")
        [g_current,fe] = grad_finitediff(fun,x_current,f_current,h);
        func_eval  = func_eval + fe;
        FE_step = FE_step + 1;
    
    
         if (FE_step == 1)
              x_pas = x_current;
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
            if ~bfgs_use
                c = 0;
            end
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

%         H =1/ (pi/20 *(cos(pi * x_current(1)/12)^2 - sin(pi * x_current(2)/16)^2));

         [beta,fval,f_temp, xi, d] =  armijo(f_current,fun,x_current,H,g_current,xi,alpha, ...
            presbounds, preslineq, presleq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas, Ai, ni, li, ui, Z, W, 1);
        
        func_eval = func_eval + fval;
        if (beta == 0 || norm(d) == 0)
            success = 0;
            evl = "LE";
            notchanged = 0;
        else
            xtemp     = x_current + beta * d/norm(d);
            g_pas     = g_current;
            x_pas     = x_current;
            y_pas     = y_current;
            success  = 1;
            x_current = xtemp;
            f_current = f_temp;
            y_current =  y_current + Z'*beta*d/norm(d);
            evl       = "FE";
            num_success_iter_fe = num_success_iter_fe + 1;


        end
        
         fvalues   = [fvalues, f_current];
         nfeval    = [nfeval, func_eval];
         tf        = [tf, fun(x_current, 1)];
         list_succ = [list_succ, success];
         list_type = [list_type,1];
    end


%     fprintf( '%5d %s %+.4e %+.4e %5d %+.4e %+.4e %.2e  %.2e\n', numIter, evl, func_eval, f_current, success, norm(g_current), alpha, xi, beta );

    count_dir = 1;%    
	if isempty(D)
		cardD = 0;
% 		warning('Iteration %d: Empty Polling set');
	else
		cardD = size(D,2);
	end
% 	numDirgen = numDirgen+cardD;

    success = 0;

    while (~success) && (count_dir <= cardD) && (func_eval < maxevals) && (evl == "LE") && (notchanged)
        %    %+.4e
		d = D(:,count_dir);
		xtemp = x_current + alpha*d;
        rho = min(1e-5, 1e-5*alpha^2);
%
%		Test for the feasibility of the point
%
		if presbounds
			xl = (min(xtemp-lb) > -tol_feas);
			xu = (max(xtemp-ub) < tol_feas);
		else
			xl = 1;
			xu = 1;
		end
%
		if ~isempty(Aeq)
			Ax = (norm(Aeq*xtemp-beq,Inf) < tol_feas);
		else
			Ax = 1;
		end
%
		if preslineq
			xli = (min(Aineq*xtemp-lbineq) > -tol_feas);
			xui = (max(Aineq*xtemp-ubineq) < tol_feas);
			Aix = (xli && xui);
		else
			Aix = 1;
		end
% %
		if (~(xl && xu && Ax && Aix))
			count_dir = count_dir +1;
		else
%
			numFeaspts = numFeaspts +1;
%
%			Test for sufficient decrease at xtemp
%
			ftemp = feval(fun,xtemp);
			func_eval = func_eval + 1;
			if (ftemp < f_current - rho*(norm(d)^2))
				success = 1;
                num_success_iter_le = num_success_iter_le + 1;

			else
				count_dir = count_dir +1;
			end
%
		end

    end
    if (evl == "LE")&& (notchanged)
        numIter = numIter + 1;
        if success
            f_current = ftemp;
            x_current = xtemp;
            y_current = y_current + alpha*Z'*d;
            alpha = gamma * alpha;
%             xi = max(alpha, xi);
%             xi = gamma * xi;
        else
            alpha = theta*alpha;
        end
        fvalues   = [fvalues, f_current];
        tf        = [tf, fun(x_current, 1)];

        nfeval    = [nfeval, func_eval];
        list_succ = [list_succ, -success];
        list_type = [list_type,-1];

        nLE       = nLE + 1;
        if (nLE > max(2,fval - 1))
            nLE = 0;
            evl = "FE";
        end
%             fprintf( '%5d %s %+.4e %+.4e  %5d  %+.4e %+.4e %.2e  %.2e\n', numIter, evl, func_eval, f_current, success, norm(g_current), alpha, xi, beta );

    end
%     end

    
    cont = (  (alpha >= tol_alpha) && (func_eval < maxevals) && (numIter < 2000*n) );
    list_beta = [list_beta, beta];






end


output.numIter = numIter;
output.funcCount = func_eval;
output.fvalues = fvalues;
output.funcCounts = nfeval;
output.FEper = FE_step;
output.truefunc = tf;
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
disp('number successful Full-Eval')
disp(num_success_iter_fe)
disp('number successful Low-Eval')
disp(num_success_iter_le)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% pause
% close all
% semilogy(list_xi, 'r')
% pause
% hold on
% semilogy(list_alpha, 'k')
% pause
