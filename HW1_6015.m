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

T = 2*pi*sqrt(42164^3/mu);
tspan = [0 T];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atApoapsisEvent);

% part a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% part b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orbitEls_f_3b, xf_3b] = hohmanTransfer(x0, NO0, dV1, dV2, 1);
sensitivity = orbitEls_f_3b-orbitEls_f_3a; % change in orbit elements per degree error

% part c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predicted_0p1 = orbitEls_f_3a + sensitivity*0.1;
predicted_0p3 = orbitEls_f_3a - sensitivity*0.3;
[orbitEls_f_0p1deg, posVel_f_0p1deg] = hohmanTransfer(x0, NO0, dV1, dV2, 0.1);
[orbitEls_f_0p3deg, posVel_f_0p3deg] = hohmanTransfer(x0, NO0, dV1, dV2, -0.3);
% Answers for part c:
diff0p1 = orbitEls_f_0p1deg - predicted_0p1;
diff0p3 = orbitEls_f_0p3deg - predicted_0p3;
% Extra Investigation:
pointErrors = linspace(-1,1,100);
error = zeros(6,length(pointErrors));
for i = 1:length(pointErrors) % get error over -1:1 deg error to plot
    [orbitEls_fi, ~] = hohmanTransfer(x0, NO0, dV1, dV2, pointErrors(i)); % simulate truth
    error(:,i) = orbitEls_fi - orbitEls_f_3a;
end
ylabels = ["a error", "e error"];
figure
for i=1:2
    subplot(2,1,i)
    plot(pointErrors, error(i,:), 'LineWidth',2)
    ylabel(ylabels(i))
    xlabel("Pointing Error (deg)")
end
sgtitle("Actual Error in Final Orbit Due to Pointing Error")


% part d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_a = 30e-3; % [km/s^2]
burnTime1 = dV1/max_a;
burnTime2 = dV2/max_a;

% orbit 1
T1 = 2*pi*sqrt(a1^3/mu);
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atPeriapsisEvent);
[t_3d1,x_3d1] = ode45(@(t,x) rdot_rddot(x,mu), tspan, x0, options); % just make the tspan 1 period maybe
xt0_3d1 = x_3d1(end,:)';
[~, NO_t0_3d] = rv2orbitEls(xt0_3d1, mu);
% transfer orbit
% dV1des_N = NO_t0_3d*[0;dV1;0];
dV1direction = NO_t0_3d*[0;1;0];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_burn1,x_burn1] = ode45(@(t,x) rdot_rddot_continuousBurn(t, x, mu, max_a, dV1direction), [t_3d1(end) t_3d1(end)+burnTime1], xt0_3d1, options);
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atApoapsisEvent);
[t_3dt,x_3dt] = ode45(@(t,x) rdot_rddot(x,mu), [t_burn1(end) t_burn1(end)+T], x_burn1(end,:)', options);
xtf_3dt = x_3dt(end,:)';
[~, NO_tf_3d] = rv2orbitEls(xtf_3dt, mu); % check that orbit elements are as expected
% orbit 2
T2 = 2*pi*sqrt(a2^3/mu);
% dV2des_N = NO_tf_3d*[0;dV2;0];
dV2direction = NO_tf_3d*[0;1;0];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_burn2,x_burn2] = ode45(@(t,x) rdot_rddot_continuousBurn(t, x, mu, max_a, dV2direction), [t_3dt(end) t_3dt(end)+burnTime2], xtf_3dt, options);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@atPeriapsisEvent);
[t_3d2,x_3d2] = ode45(@(t,x) rdot_rddot(x,mu), [t_burn2(end) t_burn2(end)+T2], x_burn2(end,:)', options);
xtf_3d2 = x_3d2(end,:)';
[orbitEls_f_3d, ~] = rv2orbitEls(x_burn2(end,:)', mu); % check that orbit elements are as expected

figure
plot3(x_3d1(:,1), x_3d1(:,2), x_3d1(:,3), 'b','LineWidth', 1.5)
hold on
plot3(x_3dt(:,1),x_3dt(:,2),x_3dt(:,3),'c', 'LineWidth', 1.5)
plot3(x_3d2(:,1),x_3d2(:,2),x_3d2(:,3),'g', 'LineWidth', 1.5)
plot3(x_burn1(:,1),x_burn1(:,2),x_burn1(:,3), 'r*','LineWidth',1)
plot3(x_burn2(:,1),x_burn2(:,2),x_burn2(:,3), 'r*','LineWidth',1)
legend("o1","o2","o3","burn 1","burn 2", "Location",'best')
grid on
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")
title("Finite Burn Trajectory (ECI)")

figure
for i=1:6
    subplot(6,1,i)
    plot(t_3d1, x_3d1(:,i))
end

finiteBurnErr = orbitEls_f_3d - orbitEls_f_3a;

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
    [r_N, v_N, r] = unpackPosVel(x);
    value = dot(r_N, v_N);
    isterminal = 1;
    % mu = 3.986e5; % [km^3/s^2] 
    % if t > 2*pi*sqrt(r^3/mu)
    %     isterminal = 1;
    % else
    %     isterminal = 0;
    % end
    direction = -1;
end

function [value,isterminal,direction] = atPeriapsisEvent(t,x)
    [r_N, v_N, r] = unpackPosVel(x);
    value = dot(r_N, v_N);
    mu = 3.986e5; % [km^3/s^2] 
    if t > 1/2*pi*sqrt(r^3/mu)
        isterminal = 1;
    else
        isterminal = 0;
    end
    direction = 1;
end

function [orbitEls_f, posVel_f] = hohmanTransfer(x0, NO0, dV1, dV2, pointingError)
% Simulate Hohman Transfer with full nonlinear equations, impulsive dVs,
% and a 2D pointing error
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


% function for ode45 to integrate
function dXdt = rdot_rddot_continuousBurn(t, X, mu, max_a, dVdirection)
        
    r = X(1:3); % extract position from state vector
    v = X(4:6); % extract velocity from state vector
    dr = v; % inertial time derivative of position is simply velocity
    dv = -mu/(norm(r)^3)*r; % inertial time derivative of velocity from equation 2.22
    
    % max_a = 30e-3; % [km/s^2]
    % dVdes = norm(dVdes_N);
    % dVdirection = dVdes_N/dVdes;
    % burnTime = dVdes/max_a;
    % if t <= burnTime
    %     dXdt = [dr; dv + max_a*dVdirection]; % assign back to state vector derivative
    % else
    %     dXdt = [dr; dv]; % assign back to state vector derivative
    % end
    
    dXdt = [dr; dv + max_a*dVdirection]; % assign back to state vector derivative

end