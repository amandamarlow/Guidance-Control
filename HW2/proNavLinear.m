function [xdot] = proNavLinear(t, x, N, tf, aug)
%PRONAVLINEAR Summary of this function goes here
%   Detailed explanation goes here

% y = x(1);
% ydot = x(2);
% nt = x(3);

A = [0, 1, 0; ...
    0, 0, 1; ...
    0, 0, 0];
B = [0; -1; 0]; 

tgo = tf-t;

nc = get_nc(x, tgo, N, aug);

xdot = A*x + B*nc;

end

