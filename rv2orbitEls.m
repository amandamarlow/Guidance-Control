function [orbitEls, NO] = rv2orbitEls(posVel, mu)
%RV2ORBITELS Summary of this function goes here
%   Detailed explanation goes here

X_N = [1;0;0];
Y_N = [0;1;0];
Z_N = [0;0;1];

r_N = posVel(1:3);
v_N = posVel(4:6);
r = norm(r_N);

h_N = cross(r_N, v_N);
h = norm(h_N);
p = h^2/mu;
e_N = cross(v_N,h_N)/mu - r_N/r;
e = norm(e_N);
a = p/(1+e^2);

% compute orbit elements
n_N = cross(Z_N, h_N);
n = norm(n_N);
OMEGA = acos(dot(n_N,X_N)/n);
if dot(n_N,Y_N) < 0
     OMEGA = -OMEGA;
end
omega = acos(dot(n_N,e_N)/n/e);
if dot(e_N,Z_N) < 0
    omega = -omega;
end
i = acos(dot(h_N,Z_N)/h);
nu = acos(dot(e_N,r_N)/e/r);
if dot(r_N,v_N) <  0
   nu = -nu; 
end
nu = real(nu);
M = atan2(-sqrt(1-e^2)*sin(nu),-e-cos(nu)) + pi-e*sqrt(1-e^2)*sin(nu)/(1+e*cos(nu));

orbitEls = [a;e;i;OMEGA;omega;M];

theta = nu + omega;
ON = R3(theta)*R1(i)*R3(OMEGA);
NO = ON';
end

