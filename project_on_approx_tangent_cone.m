function [projg, gp] = project_on_approx_tangent_cone(x_current, g_current, H, xi, preslineq, presbounds, presleq, Aeq, Ai, ni, li, ui, Z, W)
[Iap,Iam,nbgen] = activeconsAiZ(x_current,xi,li,ui,Ai,Z);
Iip = setdiff(1:ni,Iap);
Iim = setdiff(1:ni,Iam);
Ya  = [];
if nbgen == 0

%      if ~preslineq %&& ~presbounds 
        gp =  -Z*H*Z' * g_current;
        projg = -Z*Z'*g_current;
%     else
%         % Compute a direction such that Ax = b
%         gp = -Z*H*Z'*Ai'* g_current;
%         projg = -Z*Z'*Ai'*g_current;
%      end

else
    if (preslineq || (presleq && presbounds))
        [Ys,Yc,~] = gentech(Aeq,Ai,Z,Iap,Iam);
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
        ZA = Z'*Ai';
        Is = intersect(Iip,Iim);
		Icp = setdiff(Iip,Is);
		Icm = setdiff(Iim,Is);
		Ys = ZA(:,Is);
		Yc = [ZA(:,Icp) -ZA(:,Icm)];


% 		Ys = W(:,Is);
% 		Yc = [W(:,Icp) -W(:,Icm)];
    end
    dyy = 2*size(Ys,2)+size(Yc,2);    

    if (dyy > 0)
        proj_on = [Ys, -Ys, Yc];
        gp = Z*projected_gradient(g_current, proj_on, Z, H);
        %projection of minus the gradient
        projg = Z*projected_gradient(g_current, proj_on, Z, eye(size(H,1)));
    else
        n = length(x_current);
        gp = zeros(n, 1);
        projg = zeros(n, 1);
    end
        
end