% HW 3 - ASEN 6015
% Amanda Marlow
% 10/18/24

clear
clc
close all

N = 4; % length of state vector
p = 2; % length of control vector
ylabels = ["$x$ ($m$)", "$y$ ($m$)", "$\dot{x}$ ($\frac{m}{s}$)", "$\dot{y}$ ($\frac{m}{s}$)"];

mu = 3.986004418e14; % m^3/s^2
a = 6728e3; % m
n = sqrt(mu/a^3);
A = [0, 0, 1, 0;
    0, 0, 0, 1;
    3*n^2, 0, 0, 2*n;
    0, 0, -2*n, 0];
B = [0, 0;
    0, 0;
    1, 0;
    0, 1];

%% Section 1 - Natural Dynamics

% tspan = [0, 1700];
% tspan = linspace(0, 1700, 500);
tspan = linspace(0, 10000, 1000);
tf = 20*60; % s
% tspan = linspace(0, tf, 500);
X0 = [-1e3, -0.25e3, 3, -0.5]'; % m & m/s
[tnatural,Xnatural] = ode45(@(t,X) CWHode(t, X), tspan, X0);
distance = vecnorm(Xnatural(:,1:2), 2,2);
[closest,closestIdx] = min(distance);
tclosest = tnatural(closestIdx);

figure
plot(Xnatural(:,1), Xnatural(:,2), 'LineWidth',1.5)
% plot(Xnatural(:,2), Xnatural(:,1), 'LineWidth',1.5)
hold on
scatter(0,0,'filled')
scatter(Xnatural(closestIdx,1), Xnatural(closestIdx,2), 'filled')
legend("Trajectory", "ISS", "Closest Approach", 'Location','best')
axis equal
title("Natural Dynamics Trajectory")
xlabel("x (m)")
ylabel("y (m)")
% xlabel("along track (m)")
% ylabel("radial (m)")
grid on

figure
sgtitle("Natural Dynamics")
for i = 1:N
    subplot(N,1,i)
    plot(tnatural, Xnatural(:,i), 'LineWidth',1.5)
    hold on
    xline(tclosest, '--','LineWidth',1.5)
    legend('', 'time of closest approach', 'Location','northeast')
    xlabel("Time (s)", 'FontSize',12)
    ylabel(ylabels(i), 'Interpreter','latex', 'FontSize',14)
    % grid on
end

%% Section 2 - Optimal Trajectory

% Numerical with STM
A_aug = [A, -0.5*B*B'; zeros(4), -A'];
tf = 20*60; % s
Xf = [-0.25e3, 0, 0.2, 0]'; % m & m/s
Phi = expm(A_aug*tf);
phi11 = Phi(1:N,1:N);
phi12 = Phi(1:N, N+1:end);
costate0 = phi12\(Xf - phi11*X0);

time = linspace(0,tf,1000);
X_aug = zeros(2*N, length(time));
u = zeros(p, length(time));
X_aug(:,1) = [X0; costate0];
u(:,1) = -0.5*B'*X_aug(N+1:end,1);
for i = 2:length(time)
    t = time(i);
    Phi = expm(A_aug*t);
    X_aug(:,i) = Phi*X_aug(:,1);
    u(:,i) = -0.5*B'*X_aug(N+1:end,i);
end
Xoptimal = X_aug(1:4,:);
costate = X_aug(5:end,:);

figure
plot(Xoptimal(1,:), Xoptimal(2,:), 'LineWidth',1.5)
hold on
scatter(X0(1), X0(2), 'filled')
scatter(Xf(1), Xf(2), 'filled')
legend("Trajectory", "x0", "xf", 'Location','best')
axis equal
title("Optimal Trajectory")
xlabel("x (m)")
ylabel("y (m)")
grid on

figure
sgtitle("Optimal Dynamics")
for i = 1:N
    subplot(N,1,i)
    plot(time, Xoptimal(i,:), 'LineWidth',1.5)
    hold on
    scatter(time(1), X0(i), 'filled')
    scatter(time(end), Xf(i), 'filled')
    legend('', 'X0','Xf', 'Location','northeast')
    xlabel("Time (s)", 'FontSize',12)
    ylabel(ylabels(i), 'Interpreter','latex', 'FontSize',14)
    % grid on
end
figure
sgtitle("Optimal Controls")
uLabels = ["ux", "uy"];
for i = 1:p
    subplot(p,1,i)
    plot(time, u(i,:), 'LineWidth',1.5)
    xlabel("Time (s)", 'FontSize',12)
    ylabel(uLabels(i))
    grid on
end

% Part 4
[tafter,Xafter] = ode45(@(t,X) CWHode(t, X), [tf, tf+15*60], Xf);

figure
plot(Xoptimal(1,:), Xoptimal(2,:), 'LineWidth',1.5)
hold on
plot(Xafter(:,1), Xafter(:,2), 'LineWidth',1.5)
scatter(X0(1), X0(2), 'filled')
scatter(Xf(1), Xf(2), 'filled')
legend("Trajectory", "no further control", "x0", "xf", 'Location','best')
axis equal
title("Optimal Trajectory + Natural")
xlabel("x (m)")
ylabel("y (m)")
grid on
figure
sgtitle("Optimal Dynamics + Natural")
for i = 1:N
    subplot(N,1,i)
    plot(time, Xoptimal(i,:), 'LineWidth',1.5)
    hold on
    plot(tafter, Xafter(:,i), 'LineWidth',1.5)
    scatter(time(1), X0(i), 'filled')
    scatter(time(end), Xf(i), 'filled')
    % legend('optimal', 'no further control', 'X0','Xf', 'Location','northeast')
    legend('', '', 'X0','Xf', 'Location','northeast')
    xlabel("Time (s)", 'FontSize',12)
    ylabel(ylabels(i), 'Interpreter','latex', 'FontSize',14)
    % grid on
end

%% Section 3 Guidance

X0_perturbed = [-0.98e3; -0.4e3; 3.7; 0.3];
X_aug_perturbed = zeros(2*N, length(time));
u_perturbed = zeros(p, length(time));
X_aug_perturbed(:,1) = [X0_perturbed; costate0];
u_perturbed(:,1) = -0.5*B'*X_aug_perturbed(N+1:end,1);
for i = 2:length(time)
    t = time(i);
    Phi = expm(A_aug*t);
    X_aug_perturbed(:,i) = Phi*X_aug_perturbed(:,1);
    u_perturbed(:,i) = -0.5*B'*X_aug_perturbed(N+1:end,i);
end
X_perturbed = X_aug_perturbed(1:4,:);
costate_perturbed = X_aug_perturbed(5:end,:);

figure
plot(Xoptimal(1,:), Xoptimal(2,:), 'LineWidth',1.5)
hold on
plot(X_perturbed(1,:), X_perturbed(2,:), 'LineWidth',1.5)
scatter(X0(1), X0(2), 'filled')
scatter(Xf(1), Xf(2), 'filled')
legend("Optimal", "Perturbed", "x0", "xf", 'Location','best')
axis equal
title("Perturbed Trajectory")
xlabel("x (m)")
ylabel("y (m)")
grid on
figure
sgtitle("Perturbed Dynamics")
for i = 1:N
    subplot(N,1,i)
    plot(time, Xoptimal(i,:), 'LineWidth',1.5)
    hold on
    plot(time, X_perturbed(i,:), 'LineWidth',1.5)
    scatter(time(1), X0(i), 'filled')
    scatter(time(end), Xf(i), 'filled')
    legend('Optimal', 'Perturbed', 'X0', 'Xf', 'Location','northeast')
    xlabel("Time (s)", 'FontSize',12)
    ylabel(ylabels(i), 'Interpreter','latex', 'FontSize',14)
    % grid on
end


dX0 = X0_perturbed - X0;
dcostate0 = -phi12\phi11*dX0;
dX_aug = zeros(2*N, length(time));
du = zeros(p, length(time));
dX_aug(:,1) = [dX0; dcostate0];
du(:,1) = -0.5*B'*dX_aug(N+1:end,1);
for i = 2:length(time)
    t = time(i);
    Phi = expm(A_aug*t);
    dX_aug(:,i) = Phi*dX_aug(:,1);
    du(:,i) = -0.5*B'*dX_aug(N+1:end,i);
end
X_optimal_perturbed = Xoptimal + dX_aug(1:N, :);
u_optimal_perturbed = u + du;

% X_aug_neighboringOpt = zeros(2*N, length(time));
% u_neighboringOpt = zeros(p, length(time));
% X_aug_neighboringOpt(:,1) = [X0_perturbed; costate0+dcostate0];
% u_neighboringOpt(:,1) = -0.5*B'*X_aug_neighboringOpt(N+1:end,1);
% for i = 2:length(time)
%     t = time(i);
%     Phi = expm(A_aug*t);
%     X_aug_neighboringOpt(:,i) = Phi*X_aug_neighboringOpt(:,1);
%     u_neighboringOpt(:,i) = -0.5*B'*X_aug_neighboringOpt(N+1:end,i);
% end
% X_neighboringOpt = X_aug_neighboringOpt(1:4,:);

figure
plot(X_perturbed(1,:), X_perturbed(2,:), 'LineWidth',1.5)
hold on
plot(X_optimal_perturbed(1,:), X_optimal_perturbed(2,:), 'LineWidth',1.5)
% plot(X_neighboringOpt(1,:), X_neighboringOpt(2,:), 'LineWidth',1.5)
scatter(X0_perturbed(1), X0_perturbed(2), 'filled')
scatter(Xf(1), Xf(2), 'filled')
legend("Perturbed", "Neighboring Optimal", "X0 Perturbed","xf", 'Location','best')
axis equal
title("Neighboring Optimal Trajectory")
xlabel("x (m)")
ylabel("y (m)")
grid on
figure
sgtitle("Neighboring Optimal Dynamics")
for i = 1:N
    subplot(N,1,i)
    % plot(time, Xoptimal(i,:), 'LineWidth',1.5)
    hold on
    plot(time, X_perturbed(i,:), 'LineWidth',1.5)
    plot(time, X_optimal_perturbed(i,:), 'LineWidth',1.5)
    scatter(time(1), X0(i), 'filled')
    scatter(time(end), Xf(i), 'filled')
    legend('Perturbed', 'Neighboring Optimal','X0', 'Xf', 'Location','northeast')
    xlabel("Time (s)", 'FontSize',12)
    ylabel(ylabels(i), 'Interpreter','latex', 'FontSize',14)
    % grid on
end

figure
sgtitle("Neighboring Optimal Controls")
uLabels = ["ux", "uy"];
for i = 1:p
    subplot(p,1,i)
    plot(time, u(i,:), 'LineWidth',1.5)
    hold on
    plot(time, u_optimal_perturbed(i,:), 'LineWidth',1.5)
    xlabel("Time (s)", 'FontSize',12)
    ylabel(uLabels(i))
    grid on
end
legend('Nominal', 'Neighboring Optimal')