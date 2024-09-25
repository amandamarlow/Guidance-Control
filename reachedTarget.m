function [value,isterminal,direction] = reachedTarget(r, rf, tol)
%REACHEDTARGET Summary of this function goes here
%   Detailed explanation goes here
dif_r = norm(r-rf);
value = (dif_r<tol);
if value == 1
    fprintf("Reached Target")
end
isterminal = 1;
direction = 0;
end

