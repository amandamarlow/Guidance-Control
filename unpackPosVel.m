function [r_N, v_N, r, v] = unpackPosVel(x)
%UNPACKPOSVEL Summary of this function goes here
%   Detailed explanation goes here

r_N = x(1:3);
v_N = x(4:6);
r = norm(r_N);
v = norm(v_N);
end

