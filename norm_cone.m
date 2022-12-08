function [fx] = norm_cone( x,Ys, Yc, y)

gen = [Ys, -Ys, Yc];

fx = norm (gen * x - y);





