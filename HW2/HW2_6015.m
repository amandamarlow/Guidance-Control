% Homework 2 - ASEN 5090
% Amanda Marlow
% 9/24/24

clear
clc
close all

addpath("C:\Users\marlo\MATLAB Drive\6015")

%% Problem 1

Vm = 3000; % [ft/s]
Vt = 1000; % [ft/s]
R0 = 40000; % [ft]
% [t,x] = simMissile(R0, Vm, Vt, N, HE0, nt, pltTitle)
% HE0 = 30; % [deg]
% nt = 0;
% N = 3;
HE0 = [30, 0, 30];
nt = [0, 3, 3]*32; % ft/s^2
% for N = 3:5
%     for i = 1:3
%         pltTitle = strcat("N'=", num2str(N), " $\theta_{HE}$=", num2str(HE0(i)), "$^{\circ}$ $n_T$=", num2str(nt(i)), "$\frac{ft}{s^2}$");
%         [t,x] = simMissile(R0, Vm, Vt, N, HE0(i), nt(i), pltTitle);
%     end
% end

% aug = 0;
% for i = 1:3
%     figure
%     for N = 3:5
%         [t,x,nc] = simMissile(R0, Vm, Vt, N, HE0(i), nt(i), aug);
%     end
%     sgtitle(strcat("Pro Nav: ","$\theta_{HE}$ = ", num2str(HE0(i)), "$^{\circ}$ $n_T$ = ", num2str(nt(i)), "$\frac{ft}{s^2}$"), "Interpreter", "Latex")
% end
% aug = 1;
% for i = 1:3
%     figure
%     for N = 3:5
%         [t,x,nc] = simMissile(R0, Vm, Vt, N, HE0(i), nt(i), aug);
%     end
%     sgtitle(strcat("Augmented Pro Nav: ","$\theta_{HE}$ = ", num2str(HE0(i)), "$^{\circ}$ $n_T$ = ", num2str(nt(i)), "$\frac{ft}{s^2}$"), "Interpreter", "Latex")
% end

aug = [0,1];
controlEffortProNav = zeros(3,3);
controlEffortAugProNav = zeros(3,3);
missDistProNav = zeros(3,3);
missDistAugProNav = zeros(3,3);
for i = 1:3
    figure
    for j = 1:2
        for N = 3:5
            [t,x,nc] = simMissile(R0, Vm, Vt, N, HE0(i), nt(i), aug(j));
            if aug(j) == 1
                controlEffortAugProNav(N-2,i) = trapz(t,nc);
                missDistAugProNav(N-2,i) = x(end,1);
            else
                controlEffortProNav(N-2,i) = trapz(t,nc);
                missDistProNav(N-2,i) = x(end,1);
            end
        end
    end
    sgtitle(strcat("Pro Nav: ","$\theta_{HE}$ = ", num2str(HE0(i)), "$^{\circ}$ $n_T$ = ", num2str(nt(i)), "$\frac{ft}{s^2}$"), "Interpreter", "Latex")
end

% Part b
tgo_error = mvnrnd(0,0.3,100);
missDist = zeros(length(tgo_error),1);
figure
for i = 1:length(tgo_error)
    [t,x] = simMissile(R0, Vm, Vt, 4, 30, 3*32, 1, tgo_error(i));
    missDist(i) = x(end,1);
end
[t,x] = simMissile(R0, Vm, Vt, 4, 30, 3*32, 1);

figure
scatter(tgo_error, missDist)
yline(0,'k--')
xlabel("Time To Go Error (s)")
ylabel("Miss Distance (m)")
title("Monte Carlo Runs with Error in Time to Go")

%% Problem 2 - Lambert guidance
r0 = 6578;
r0_N = [r0;0;0];
rf = 42164;
rf_N = [-rf*cosd(5); rf*sind(5); 0];
mu = 3.986e5;

TOF = 37864;
using_long_transfer_angle = 0;
[v0_N_a, vf_N, iterations] = Lamberts(mu, r0_N, rf_N, using_long_transfer_angle, TOF);

% Part B
v0 = sqrt(mu/r0);
v0_N = [0;v0;0];
x0_N = [r0_N; v0_N];
at = 30/1000; % km/s^2
[t_hist, x_hist, at_hist] = LambertGuidance(x0_N, rf_N, TOF, at, mu);
error_N = x_hist(end,1:3)'-rf_N;
error = norm(error_N);

firstBurn = trapz(t_hist(1:140),double(at_hist(1:140)>0));
burnTime_tot = trapz(t_hist,double(at_hist>0));

% plotting
figure
scatter(x_hist(:,1),x_hist(:,2), '.')
hold on
scatter(rf_N(1),rf_N(2),'*', 'LineWidth',4)
axis equal
grid on
xlabel("x [km]")
ylabel("y [km]")
title("Lambert Guidance Trajectory")

% figure
% plot(t_hist,at_hist)
% ylim([-0.001, 0.031])
figure
scatter(t_hist,at_hist, 'filled')
ylim([-0.001, 0.031])
ylabel("acceleration [km/s^2]")
xlabel("time [s]")
title("Acceleration vs. Time - Lambert Guidance")

% figure
% plot3(x_hist(:,1),x_hist(:,2),x_hist(:,3), '.')
% hold on
% scatter3(rf_N(1),rf_N(2),rf_N(3),'*', 'LineWidth',4)
% axis equal
% grid on

% figure
% for i=1:6
%     subplot(6,1,i)
%     plot(t_hist, x_hist(:,i))
% end

% Part C
crossProdSteer = 1;
[t_hist_c, x_hist_c, at_hist_c] = LambertGuidance(x0_N, rf_N, TOF, at, mu, crossProdSteer);
error_N_C = x_hist_c(end,1:3)'-rf_N;
error_C = norm(error_N_C);
firstBurn_C = trapz(t_hist_c(1:140),double(at_hist_c(1:140)>0));
burnTime_tot_C = trapz(t_hist_c,double(at_hist_c>0));

% plotting
figure
scatter(x_hist_c(:,1),x_hist_c(:,2), '.')
hold on
scatter(rf_N(1),rf_N(2),'*', 'LineWidth',4)
axis equal
grid on
xlabel("x [km]")
ylabel("y [km]")
title("Cross Product Steering Trajectory")

figure
scatter(t_hist,at_hist)
hold on
scatter(t_hist_c,at_hist_c, '.')
legend("Lambert", "Cross Product")
ylabel("acceleration [km/s^2]")
xlabel("time")
title("Acceleration vs. Time Comparison")