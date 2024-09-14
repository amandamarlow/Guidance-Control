function [NO, ON] = DCMfromOrbitEls(orbitEls)
%DCMFROMORBITELS Summary of this function goes here
%   Detailed explanation goes here
[~, e, i, OMEGA, omega, ~] = unpackOrbitEls(orbitEls);

% nu = acos(dot(e_N,r_N)/e/r);
theta = nu + omega;
ON = R3(theta)*R1(i)*R3(OMEGA);
NO = ON';
end

