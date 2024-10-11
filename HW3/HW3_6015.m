% HW 3 - ASEN 6015
% Amanda Marlow
% 10/8/24

clear
clc
close all

%% Question 1

A = [5 0 3 0 1;
    3 0 0 -2 0;
    0 -2 4 1 0;
    1 3 -4 1 3;
    0 2 2 0 -1];
B = [0 1;
    0 2;
    0 0;
    1 3;
    1 1];
C = [0 0 1 0 0];

gamma = 1;
a = 1;
b = 1;
Q = zeros(5,5);
Q(3,3) = gamma;
R = [a, 0; 0, b];
[K,S,P] = lqr(A,B,Q,R);
sys_LQR = ss(A-B*K,B,C, zeros(1,2));
IC = ones(5,1);
[y_LQR,t_LQR,x_LQR] = initial(sys_LQR,IC);
figure
plot(t_LQR,y_LQR, 'LineWidth',1.5)
title("Q1: Output Response")
figure
for i = 1:5
 subplot(5,1,i)
 plot(t_LQR, x_LQR(:,i), 'LineWidth',1.5)
end
sgtitle("Q1: State Response")

%% Question 2

H = A - B*K;
R_wi = [0.1, 0.2, 0.1, 0.3, 0.2];
Rw = diag(R_wi);
sys_adjoint = ss(-(A-B*K),zeros(5,1),zeros(1,5), zeros(1,1));
p0 = C';
% [p,tOut,x] = initial(sys_adjoint,IC,flip(t_LQR));
tspan = [t_LQR(end), 0];
[t,p] = ode45(@(t,x) adjointODE(t, x, H), tspan, p0);
figure
plot(t,p, 'LineWidth',1.5)

integrand = zeros(size(t));
integrand_channels = zeros(length(t),5);
for i = 1:length(t)
    integrand(i) = p(i,:)*Rw*p(i,:)';
    for j = 1:5
        integrand_channels(i,j) = R_wi(j)*p(i,j)^2;
    end
end
var_y = trapz(flip(t),flip(integrand));
var_channels = trapz(flip(t),flip(integrand_channels,1),1);

% part B
% gamma = 200;
% Q = zeros(5,5);
% Q(3,3) = gamma;

% Kf = zeros(5,5); % S
Kf = zeros(25,1); % S
tspan = linspace(t_LQR(end),0,1000);
% tspan = linspace(10,0,1000);
[t_hist,K] = ode45(@(t,K) K_ODE(t, K, A, B, Q, R), tspan, Kf);
K = reshape(K',5,5,[]);

x0 = ones(5,1);
K_hist = permute(K,[3,1,2]);
% tspan = [0, 10];
tspan = flip(tspan);
[t,x] = ode45(@(t,x) timeVaryingLQR(t, x, A, B, R, t_hist, K_hist), tspan, x0);

figure
plot(t,x, 'LineWidth',1.5)
legend("x1", "x2", "x3", "x4", "x5")
figure
plot(t,x(:,3), 'LineWidth',1.5)

% x0 = ones(5,1);
% S0 = [x0;reshape(K0,[],1)];
% [t,S] = ode45(@(t,K) xK_ODE(t, S, A, B, Q, R), t_LQR, S0);
% x = S(:,1:5);
% K = reshape(S(:,6:end)',5,5,[]);
% STM_vec = reshape(flatSTM',n,n,[]);