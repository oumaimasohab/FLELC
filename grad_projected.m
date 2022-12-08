function [g,func_eval] = grad_projected(x,f,h,D)
n = length(size(D,2));
g = zeros(length(x),1);

for i=1:n
    D(:,i) = D(:,i)/norm(D(:,i));
    try
        [fval, ~] = spamfunc(x+h*D(:,i), 0);
    catch
        warning('grad approx has trouble computing f. exiting program.');
        pause
        return
    end 
    g = g +  max((fval-f)/h, 0) * D(:,i);   
end
func_eval = n;