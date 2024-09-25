function [dXdt] = LambertGuidance(t, x, rf_N, TOF, at, mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r_N = x(1:3);
r = norm(r_N);
v_N = x(4:6);
tgo = TOF-t;

tol = 1e-6; % km/s
using_long_transfer_angle = 0;
[vr_N, ~, ~] = Lamberts(mu, r_N, rf_N, using_long_transfer_angle, tgo);
vg_N = vr_N - v_N;
vhat = vg_N/norm(vg_N);
if all(vg_N < ones(3,1)*tol)
    at_N = zeros(3,1);
else
    at_N = at*vhat;
end

dr = v_N; % inertial time derivative of position is simply velocity
dv = -mu/(norm(r)^3)*r_N + at_N; % inertial time derivative of velocity from equation 2.22
dXdt = [dr; dv]; % assign back to state vector derivative
end

