function [f] = fun(x, istrue, varargin)

probtype = '';
try
    [f, ~] = spamfunc(x, 0);
catch
    f = Inf;
end

if (nargin<2)
    istrue = 0;
end

if (~istrue)
    switch probtype
        case 'noisy3'
            sigma=10^-3;
            u = sigma*(-1+2*rand(1,1));
            f = f*(1+u);
%             f = f + u;
        case 'wild3'
            sigma=10^-3;
            phi = 0.9*sin(100*norm(x,1))*cos(100*norm(x,inf)) + 0.1*cos(norm(x,2));
            phi = phi*(4*phi^2 -3);
            f = (1 + sigma*phi)*f;
%             f = sigma*phi + f;
    end
end