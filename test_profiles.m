list_problems = get_list_problems(["linearequ_const_copy.csv", "linearinequ_const_less.csv"]);
nsample = 1;
H = zeros(nsample,nsample,3);
N = zeros(nsample,nsample,3);
S = zeros(nsample,nsample);
T = zeros(nsample,nsample);


Ngrads = zeros(length(list_problems)*nsample);

for i=1:length(list_problems)
    prob          = list_problems(i);  
    p             = prob{:};
    disp(['problem ' int2str(i) ': ' char(p)])
    [x0, lb, ub, Aeq, beq, Aineq, lbineq, ubineq] = get_problem(p);
    n = length(x0);

    dim = 2000*n;
    tol_feas = 1e-8;

%     if (n < 5)
%         i
%     end

%     
    for s=1:nsample
        if (i==20)
            prob
            Ngrads(nsample*(i-1)+s) = length(x0)+1;
%         if ( i== 35)
            [xsol1,fsol1,~,output_ds] = dspfd(@fun,x0,1,2,lb,ub, Aeq, beq, Aineq, lbineq, ubineq);
            truef = output_ds.truefunc;
            hist = output_ds.histF;
            fvalues0 = hist(:,1);
            nnfev0 = hist(:,2);
            H(1:length(fvalues0),nsample*(i-1)+s,1) = truef;
            N(1:length(fvalues0),nsample*(i-1)+s,1) = nnfev0;

             [xsol2,fsol2,~,output_fe, typeofeval,  list_xi1, list_alpha1, list_beta1, list_successful, list_type] = FLE_const(@fun,x0,lb,ub,Aeq, beq,Aineq,lbineq,ubineq,dim,tol_feas, 1, 0);
            fvalues1 = output_fe.fvalues;
            nnfev1 = output_fe.funcCounts;
            truef1 = output_fe.truefunc;

            H(1:length(fvalues1),nsample*(i-1)+s,2)   = truef1;
            N(1:length(fvalues1),nsample*(i-1)+s,2)   = nnfev1;
            S(1:length(fvalues1)-1,nsample*(i-1)+s)   = list_successful;
            T(1:length(fvalues1)-1,nsample*(i-1)+s)   = list_type;
        end
% 
%             
%             opts = nomadset('bb_output_type',bbtype);
%             opts.max_bb_eval = 2000*length(x0);
%             opts.model_search = 1;
%             opts.stats_file = './output.txt';
%             opts.history_file = './iterates.txt';
% 
% 
%     
%             [x,fval,exitflag,iter,nfval] = nomad(@(x)fun_nomad(x, A, d, e), x0, lb, ub,opts);
%     
%             fileID = fopen('/home/ouma/Documents/rBFGS/constrained_LF/output.0.txt','r');
%             formatSpec = '%d %f';
%             sizeA = [2 Inf];
%             feval = fscanf(fileID,formatSpec,sizeA);
%             fclose(fileID);
% 
%             fileID_iterates = fopen('/home/ouma/Documents/rBFGS/constrained_LF/iterates.0.txt','r');
%             iterates = fscanf(fileID_iterates, repmat('%f ', [1, n]), [n+1 Inf]);
%             fvalues = [];
%             if (~isempty(feval))
%                 for ev=feval(1,:)
%                     x = iterates(1:n, ev);
%                     fvalues = [fvalues, fun(x, 1)];
%                 end
%                 H(1:size(feval,2),nsample*(i-1)+s,3) = fvalues;
%                 N(1:size(feval,2),nsample*(i-1)+s,3) = feval(1,:);
%             end



%         end
%         pause
% %         pause
%         [~,fsol2,~,output_fe, typeofeval,  list_xi2, list_alpha2] = Full_eval_const(@fun,x0,lb,ub,Aeq, beq,Aineq,lbineq,ubineq, dim, tol_feas); 
% %
% %           [~,fsol2,~,output_fe, typeofeval,  list_xi2, list_alpha2, list_beta2] = Low_eval_const(@(x)fun(x),x0,lb,ub,Aeq, beq,Aineq,lbineq,ubineq,dim,tol_feas, 1, 0);
%         fvalues3 = output_fe.fvalues;
%         nnfev3 = output_fe.funcCounts;
%         truef3 = output_fe.truefunc;
% 
%         H(1:length(fvalues3),nsample*(i-1)+s,4) = truef3;
%         N(1:length(fvalues3),nsample*(i-1)+s,4) = nnfev3;
        
%         close all
% 
%         figure(1)
% 
%         plot(fvalues0, 'b')
%         hold on
%         plot(fvalues1,'r')
%         plot(fvalues2, 'k')
%         pause

% 
%         lb = a;
%         ub = b;
%         opts = nomadset(); 
%         opts.max_bb_eval = 2000*length(x0);
%         opts.model_search = 1;
%         opts.stats_file = './output.txt';
% 
%         [x,fval,exitflag,iter,nfval] = nomad(@(x)fun(x),x0,lb,ub,opts);

% 
%         fileID = fopen('/home/ouma/Documents/rBFGS/constrained_LF/output.0.txt','r');
%         formatSpec = '%d %f';
%         sizeA = [2 Inf];
%         feval = fscanf(fileID,formatSpec,sizeA);
%         fclose(fileID);
%          H(1:size(feval,2),nsample*(i-1)+s,4) = feval(2,:);
%          N(1:size(feval,2),nsample*(i-1)+s,4) = feval(1,:);
%         [~,fsol2,~,output_fe, typeofeval,  list_xi2, list_alpha2] = Full_eval_const(@(x)fun(x),x0,a,b,[],[],[],[],[],dim,tol_feas, 1, 1);
%         fvalues3 = output_fe.fvalues;
%         nnfev3 = output_fe.funcCounts;
%         H(1:length(fvalues3),nsample*(i-1)+s,4) = fvalues3;
%         N(1:length(fvalues3),nsample*(i-1)+s,4) = nnfev3;
%         pause

%         semilogy(list_alpha1, '-r')
%         hold on
%         semilogy(list_xi1, '-b')
%         hold on
%         semilogy(list_beta1, '-k')
% %         xlim([0 min(100, length(fvalues1))])
%         ylim([0 max(list_xi1)])
% 
% 
%         legend('\alpha', '\xi', '\beta')
%         title("Comparision of \alpha, \xi, and \beta")
% %         pause
%         saveas(gcf,['alpha_xi_values/' p '_nobfgs_500xi.fig'])
%         close
% 
% 
%         semilogy(list_alpha2, '-r')
%         hold on
%         semilogy(list_xi2, '-b')
%         hold on
%         semilogy(list_beta2, '-k')
% 
% %         xlim([0 min(100, length(fvalues2))])
%         ylim([0 max(list_xi2)])
% 
% 
%         legend('\alpha', '\xi', '\beta')
%         title("Comparision of \alpha, \xi, and \beta")
% %         pause
%         saveas(gcf,['alpha_xi_values/' p '_bfgs_500xi.fig'])
%         close
    end  
    
    
end

Ind = find(~N);

N(Ind) = inf;
H(Ind) = inf;
figure(1)
perf_profile(H, N, S, T, 1e-3, 1)
% perf_profile(H, N, S, T, 1e-5,  1)
% 
% file = 'comp_FE_version_28Feb_5xi';
% save(file,'H','N')
