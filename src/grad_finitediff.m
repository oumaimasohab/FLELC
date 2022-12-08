function [g,func_eval] = grad_finitediff(fun, x,f,h)
n = length(x);
g = zeros(n,1);
I = eye(n);

for i=1:n
    try
        fval = feval(fun, x+h*I(:,i));
    catch
        warning('grad approx has trouble computing f. exiting program.');
        return
    end 
    g(i) = (fval-f)/h;   
end
func_eval = n;

% fval_vect = [];
% G = [];
% for i=1:n
%     try
%         fval1 = feval(fun, x+h*I(:,i));
%         fval2 = feval(fun, x-h*I(:,i));
% 
%     catch
%         warning('grad approx has trouble computing f. exiting program.');
%         return
%     end 
% 
%     fval_vect = [fval_vect; fval1; fval2];
% 
%     z = zeros(1, n);
% 
%     z(i) = 1;
% 
%     G = [G ; z];
%     z(i) = -1;
%     G = [G ; z];
% 
% end
% 
% g = G\fval_vect;
% g = g/h;
% func_eval = 2*n;

