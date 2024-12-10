clear
clc
close all
profile on
addpath("C:\Users\marlo\MATLAB Drive\6015")

%% Problem Setup
% canonical units are defined from earth's surface
TU = 806.812; % s
DU = 6378.140; % km
mu_canonical = 1; % canonical units
mu = 3.986004415e5; % km^3/s^2

Fmax = 0.01; % canonical units (9.8e-5 km/s^2)

r0 = [-0.70545852988580, -0.73885031681775, -0.40116299069586]'; % canonical units
v0 = [0.73122658145185, -0.53921753373056, -0.29277123328399]'; % canonical units
x0 = [r0; v0];
[orbitEls_t0, ~] = rv2orbitEls(x0, mu_canonical, 1);

L_targ = [0, 0, 2.56612389857378]'; % Angular momentum vector of target orbit
A_targ = [0, 0, 0]'; % Laplace vector of target orbit (mu * eccentricity)
orbitEls_targ = [42000/DU; 0; 0; 0; 0; 0];


%% Tuning Parameters
epsilon = 0.00001;
k = 2; % k determines the relative weighting between the two quadratic terms in the function V

%% Simulate Dynamics

% t_max = 20;
% t_max = 18.87/TU*60^2;
t_max = 20/TU*60^2;
tspan = linspace(0,t_max,1000);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @(t,x)reachedTarget(t, x, L_targ, A_targ, mu_canonical));
% options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @(t,x)Fsmall(t, x, Fmax, L_targ, A_targ, k, epsilon, mu_canonical));
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_canonical,x_canonical] = ode45(@(t,x) lowThrustDynamics(x, L_targ, A_targ, Fmax, k, epsilon, mu_canonical), tspan, x0, options);
x_canonical = x_canonical';
n = size(x_canonical,2);
F_canonical = NaN(3,n);
G_canonical = NaN(3,n);
for i = 1:n
    [F_canonical(:,i), G_canonical(:,i)] = control_F(x_canonical(1:3,i), x_canonical(4:6,i), L_targ, A_targ, Fmax, k, epsilon, mu_canonical);
end
% [F_canonical, G_canonical] = control_F(x_canonical(1:3,:), x_canonical(4:6,:), L_targ, A_targ, Fmax, k, epsilon, mu_canonical);
x = [x_canonical(1:3,:)*DU; x_canonical(4:6,:)*DU/TU];
t = t_canonical*TU;
F = F_canonical*DU/(TU^2);

%% Plotting and Analysis
tf = t(end);
xf = x(:, end);
[orbitEls, NO] = rv2orbitEls(xf, mu, 1);
[Lf, Af] = x2LA(xf, mu);
error = [Lf - L_targ; Af - A_targ];

fprintf("final time = %f\n", tf)
fprintf("error in [L; A]:\n")
disp(error)
fprintf("af = %f, ef = %f, if = %f \n", orbitEls(1), orbitEls(2), orbitEls(3))

figure
plot3(x_canonical(1,:), x_canonical(2,:), x_canonical(3,:), "LineWidth", 1.5)
hold on
plotOrbit(100, orbitEls_t0, "--", 1.5)
plotOrbit(100, orbitEls_targ, "--", 1.5)
scatter3(x0(1), x0(2), x0(3),'k')
scatter3(x_canonical(1,end), x_canonical(2,end), x_canonical(3,end),'k*')
legend("trajectory", "initial orbit", "final orbit", "Location","best")
xlabel("x (canonical units)")
ylabel("y (canonical units)")
zlabel("z (canonical units)")
grid on
axis equal

figure
plot(t/60^2,F(1,:)*1e3, "LineWidth",2)
hold on
plot(t/60^2,F(2,:)*1e3, "LineWidth",2)
plot(t/60^2,F(3,:)*1e3, "LineWidth",2)
legend("Fx", "Fy", "Fz")
xlabel("time (hours)")
ylabel("Control Force (m/s^2)")

profile off
profile viewer
%% Functions
function [value,isterminal,direction] = reachedTarget(t, x, L_targ, A_targ, mu)
    % r = x(1:3);
    % v = x(4:6);
    % L = cross(r, v);
    % A = cross(v, L) - mu*r/norm(r);
    [L, A] = x2LA(x, mu);
    delta_L = L - L_targ;
    delta_A = A - A_targ;
    value = norm(delta_L) + norm(delta_A);
    isterminal = 1;
    direction = -1;
end

function [value,isterminal,direction] = Fsmall(t, x, Fmax, L_targ, A_targ, k, epsilon, mu)
    r = x(1:3);
    v = x(4:6);
    % L = cross(r, v);
    % A = cross(v, L) - mu*r/norm(r);
    % [L, A] = x2LA(x, mu);
    F = control_F(r, v, L_targ, A_targ, Fmax, k, epsilon, mu);
    value = norm(F) - epsilon*Fmax;
    isterminal = 1;
    direction = -1;
end

function [L, A] = x2LA(x, mu)
    r = x(1:3);
    v = x(4:6);
    L = cross(r, v);
    A = cross(v, L) - mu*r/norm(r);
end

function [] = plotOrbit(num_points, orbitEls, lineStyle, lineWidth)
    % Preallocate arrays for Cartesian coordinates
    [a, e, i, OMEGA, omega, ~] = unpackOrbitEls(orbitEls);
    % True anomaly range
    nu_range = linspace(0, 2*pi, num_points);

    x = zeros(1, num_points);
    y = zeros(1, num_points);
    z = zeros(1, num_points);
    
    % Calculate position for each point on the orbit
    for j = 1:num_points
        % Current true anomaly
        nu_current = nu_range(j);
        
        % Radius at current true anomaly
        r = (a * (1 - e^2)) / (1 + e * cos(nu_current));
        
        % Position in the orbital plane
        x_orbital = r * cos(nu_current);
        y_orbital = r * sin(nu_current);
        
        % Convert to 3D Cartesian coordinates
        x(j) = (cos(omega) * cos(OMEGA) - sin(omega) * sin(OMEGA) * cos(i)) * x_orbital + (-sin(omega) * cos(OMEGA) - cos(omega) * sin(OMEGA) * cos(i)) * y_orbital;
        y(j) = (cos(omega) * sin(OMEGA) + sin(omega) * cos(OMEGA) * cos(i)) * x_orbital + (-sin(omega) * sin(OMEGA) + cos(omega) * cos(OMEGA) * cos(i)) * y_orbital;
        z(j) = (sin(omega) * sin(i)) * x_orbital + (cos(omega) * sin(i)) * y_orbital;
    end
    
    % Plot the orbit in 3D
    plot3(x, y, z, lineStyle, "lineWidth", lineWidth);
end