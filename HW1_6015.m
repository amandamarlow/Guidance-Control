% HW 1 - ASEN 6015
% Amanda Marlow
% 9/10/24

clear
clc
close all

addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")

% %% Problem 2
% 
% a0 = 7000; % [rad]
% e0 = 0.1; % [rad]
% i0 = .7; % [rad]
% OMEGA0 = 1; % [rad]
% omega0 = 1; % [rad]
% M0 = 1; % [rad]
% 
% x0 = [a0; e0; i0; OMEGA0; omega0; M0];
% xf = [a0; e0; 1; OMEGA0; omega0; M0];
% % dxdt = GaussVariationalEq(0, x0, xf);
% tspan = [0 43200];
% % tspan = [0 10000];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12); % given tolerance
% [t_hist,x_hist] = ode45(@(t,x) proportionalControl(t, x, xf), tspan, x0, options);
% x_hist(:,6) = wrapTo2Pi(x_hist(:,6));
% 
% for i=1:6
%     subplot(6,1,i)
%     plot(t_hist,x_hist(:,i))
% end

%% Problem 3

mu = 3.986e5;
alt1 = 200;
alt2 = 35786;
Re = 6378;
a1 = alt1+Re;
a2 = alt2+Re;

at = (a1+a2)/2;
et = 1 - a1/at;
dV1 = sqrt(mu/at*(1+et)/(1-et)) - sqrt(mu/a1);
dV2 = sqrt(mu/a2) - sqrt(mu/at*((1-et)/(1+et)));

% simulate Hohman transfer
% e0 = 1e-13;
% x0 = [a1; e0; 0.1; 0; 0; 0];
% [~, B] = GaussVariationalEqs(x0);
% xt0 = x0 + B*[0; dV1; 0];

% orbitEls_0 = [a1; 0; 0.01; 0; 0; 0];
orbitEls_0 = [a1; 0; 0.01; 0; 0; 0];
[x0, NO0] = orbitEls2xv(orbitEls_0, mu);

tspan = [0 2*pi*sqrt(42164^3/mu)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atApoapsisEvent);

% part a
xt0_3a = x0 + [zeros(3,1); NO0*[0;dV1;0]];
% [orbitEls_t0_3a] = rv2orbitEls(x0, mu); % check that orbit elements are as expected

[t_3a,x_3a] = ode45(@(t,x) rdot_rddot(x,mu), tspan, xt0_3a, options);
xtf_3a = x_3a(end,:)';
[orbitEls_tf_3a, NOf_3a] = rv2orbitEls(xtf_3a, mu); % check that orbit elements are as expected
% NOf_3a = DCMfromOrbitEls(orbitEls_tf_3a);
xf_3a = xtf_3a + [zeros(3,1); NOf_3a*[0;dV2;0]];
[orbitEls_f_3a] = rv2orbitEls(xf_3a, mu); % check that orbit elements are as expected
% [~, B] = GaussVariationalEqs(xtf_3a);
% xf_3a = xtf_3a + B*[0; dV2; 0];

% part b
% xt0_3b = x0 + [zeros(3,1); NO0*[dV1*sind(1);dV1*cosd(1);0]];
% [t_3b,x_3b] = ode45(@(t,x) rdot_rddot(x,mu), tspan, xt0_3b, options);
% xtf_3b = x_3b(end,:)';
% [~, NOf_3b] = rv2orbitEls(xtf_3b, mu); % get DCM
% xf_3b = xtf_3b + [zeros(3,1); NOf_3b*[0;dV2;0]];
% [orbitEls_f_3b] = rv2orbitEls(xf_3b, mu); % check that orbit elements are as expected
[orbitEls_f_3b, xf_3b] = hohmanTransfer(x0, NO0, dV1, dV2, 1);

% part d
sensitivity = orbitEls_f_3b-orbitEls_f_3a; % change in orbit elements per degree error
predicted_0p1 = orbitEls_f_3a + sensitivity*0.1;
predicted_0p3 = orbitEls_f_3a - sensitivity*0.3;
[orbitEls_f_0p1deg, posVel_f_0p1deg] = hohmanTransfer(x0, NO0, dV1, dV2, 0.1);
[orbitEls_f_0p3deg, posVel_f_0p3deg] = hohmanTransfer(x0, NO0, dV1, dV2, -0.3);
diff0p1 = orbitEls_f_0p1deg - predicted_0p1;
diff0p3 = orbitEls_f_0p3deg - predicted_0p3;

pointErrors = linspace(-1,1,100);
error = zeros(6,length(pointErrors));
for i = 1:length(pointErrors)
    [orbitEls_fi, ~] = hohmanTransfer(x0, NO0, dV1, dV2, pointErrors(i));
    error(:,i) = orbitEls_fi - orbitEls_f_3a;
end
ylabels = ["a error", "e error"];
for i=1:2
    subplot(2,1,i)
    plot(pointErrors, error(i,:))
end

%% function
function dxdt = proportionalControl(t, x, xd)
    K = [0, 0, 0, 0,  0, 0; ...
        0, 0, 0, 0, 0, 0; ...
        0, 0, 1, 0, 0, 0];
    u = -K*(x-xd);
    
    [fx, B] = GaussVariationalEqs(x);

    dxdt = fx + B*u;
end

function dxdt = noControl(t, x)
    u = zeros(3,1);  
    [fx, B] = GaussVariationalEqs(x);
    dxdt = fx + B*u;
end

function [fx, B] = GaussVariationalEqs(x)
    mu = 3.986e5; % [km^3/s^2] 
    
    [a, e, i, OMEGA, omega, M] = unpackState(x);
    M = wrapTo2Pi(M);
    nu = nuFromM(M,a,e,mu);
    p = a*(1-e^2);
    h = sqrt(mu*p);
    r = p/(1+e*cos(nu));
    n = sqrt(mu/a^3);
    theta = omega + nu;
    
    fx = [zeros(5,1); n];
    B = [2*(a^2)/h*e*sin(nu), 2*(a^2)*p/h/r, 0; ...
        p/h*sin(nu), ((p+r)*cos(nu)+r*e)/h, 0;...
        0, 0, r*cos(theta)/h;...
        0, 0, r*sin(theta)/h/sin(i);...
        -p*cos(nu)/h/e, (p+r)/h/e*sin(nu), -r*sin(theta)*cos(i)/h/sin(i);...
        a*sqrt(1-e^2)/a/h/e*(p*cos(nu)-2*r*e), -a*sqrt(1-e^2)/a/h/e*(p+r)*sin(nu), 0];
end

% function M = M_from_nu(nu,x)
%     [~,e] = unpackState(x);
%     M = atan2(-sqrt(1-e^2)*sin(nu),-e-cos(nu)) + pi-e*sqrt(1-e^2)*sin(nu)/(1+e*cos(nu));
% end

function [nu] = nuFromM(M,a,e,mu)
    % nu = fsolve(@(nu) M-M_from_nu(nu,x), M);
    % [a,e] = unpackState(x);
    E = M2E(M, a, e, mu);
    nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end

function [a, e, i, OMEGA, omega, M] = unpackState(x)
        a = x(1);
        e = x(2);
        i = x(3);
        omega = x(4);
        OMEGA = x(5);
        M = x(6);
end

function [value,isterminal,direction] = atApoapsisEvent(t,x)
%     value = y(6)-pi; %mean anomaly at 180 degrees hopefully?
%     isterminal = 1;
%     direction = 0;
    [r_N, v_N] = unpackPosVel(x);
    value = dot(r_N, v_N);
    isterminal = 1;
    direction = -1;
end

function [orbitEls_f, posVel_f] = hohmanTransfer(x0, NO0, dV1, dV2, pointingError)
mu = 3.986e5;
xt0 = x0 + [zeros(3,1); NO0*[-dV1*sind(pointingError);dV1*cosd(pointingError);0]];
tspan = [0 2*pi*sqrt(42164^3/mu)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atApoapsisEvent);
[~,x] = ode45(@(t,x) rdot_rddot(x,mu), tspan, xt0, options);
xtf = x(end,:)';
[~, NOf] = rv2orbitEls(xtf, mu); % get DCM
posVel_f = xtf + [zeros(3,1); NOf*[0;dV2;0]];
[orbitEls_f] = rv2orbitEls(posVel_f, mu); % check that orbit elements are as expected
end