function [posVel, NO] = orbitEls2xv(orbitEls, mu)
%ORBITELS2XV Summary of this function goes here
%   NO is the DCM from orbit frame to inertial

% addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")

[a, e, i, OMEGA, omega, M] = unpackOrbitEls(orbitEls);

p = a*(1-e^2);
h = sqrt(mu*p);
nu = nuFromM(M,a,e,mu);
r = p/(1+e*cos(nu));
r_O = [r; 0; 0];
v_O = mu/h*[e*sin(nu); 1+e*cos(nu); 0];

theta = omega + nu;
% ON = Euler3132C([OMEGA;i;theta]);
ON = R3(theta)*R1(i)*R3(OMEGA);
NO = ON';
r_N = NO*r_O;
v_N = NO*v_O;
posVel = [r_N; v_N];
end



