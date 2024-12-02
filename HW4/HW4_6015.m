
clear
clc
close all
addpath("C:\Users\marlo\MATLAB Drive\6015")

% constants
mu = 3.986004415e5;

x0 = [8276; 5612; 5; -3.142; 4.672; 0];
% initial guesses
A_guess = [1; 0; 0];
B_guess = [0; 0; 0];
tf_est = 23; % burnout time
% [~, NO_t0] = rv2orbitEls(x0, mu);
% A_N = NO_t0*A_O;
% B_N = NO_t0*B_O;

a = 10000; % km
inc = 0; 
e = 0.0001; % km
omega = 35*pi/180; % deg
targetOrbitEls = [a; e; inc; 0; omega; 0];

[x_target, ~] = orbitEls2xv(targetOrbitEls, mu);
% [orbitEls0, NO_0] = rv2orbitEls(x0, mu);

rp = a*(1-e);
vp = sqrt(2*mu/rp - mu/a);
rp_O = [rp; 0; 0];
vp_O = [0; vp; 0];
[~, NO_periapsis] = orbitEls2xv([a; e; inc; 0; omega; 0], mu);
rp_N = NO_periapsis*rp_O;
vp_N = NO_periapsis*vp_O;
rp_hat = rp_N/rp;
vp_hat = vp_N/vp;
h = sqrt(mu*a*(1-e^2));
h_N = cross(rp_N, vp_N);
h_hat = h_N/h;

% p = a*(1-e^2);
p = vp^2*rp^2/mu;

maxThrust = 30e-3; % km/s



%% Question 1

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A_guess, B_guess, maxThrust), [0 tf_est], x0, options);
% T_target = 2*pi*sqrt(a^3/mu);
% [~,x_target] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A_guess, B_guess, 0), [0 T_target], x_target, options);

linewidth = 1.5;
figure
subplot(2,1,1)
plot(t,x(:,1:3), "LineWidth",linewidth)
ylabel("position (km)")
xlabel("time (s)")
legend("$r_x$","$r_y$","$r_z$", 'interpreter', 'latex')
subplot(2,1,2)
plot(t,x(:,4:6), "LineWidth",linewidth)
ylabel("velocity (km/s)")
xlabel("time (s)")
legend("$v_x$","$v_y$","$v_z$", 'interpreter', 'latex')
sgtitle("Question 1: Trajectory with no guidance")

figure
plot3(x(:,1), x(:,2), x(:,3), "r", "linewidth", 2)
hold on
% plot3(x_target(:,1), x_target(:,2), x_target(:,3), "b--","linewidth", 1)
xlabel("x (km)")
ylabel("y (km)")
zlabel("z (km)")
grid on

x_tf_N = x(end,:)';
r_tf_N = x(end,1:3)';
r_tf = norm(r_tf_N);
v_tf_N = x(end,4:6)';


%% Find best finite differencing values
lambdaAug = [B_guess; A_guess; tf_est];
t = 0;
E = ode45wrapper(t, x0, lambdaAug, maxThrust, targetOrbitEls, mu);

exponent = linspace(-12, 0, 50);
diff_test = (10*ones(size(exponent))).^exponent;
variance = NaN(size(diff_test, 2)-1, 7);
Ji = jacobian(@(input) ode45wrapper(t, x0, input, maxThrust, targetOrbitEls, mu), lambdaAug, diff_test(1)*ones(7,1), E);
for i = 1:size(diff_test, 2)-1
    Jii = jacobian(@(input) ode45wrapper(t, x0, input, maxThrust, targetOrbitEls, mu), lambdaAug, diff_test(i+1)*ones(7,1), E);
    Jdiff = Jii-Ji;
    variance(i,:) = vecnorm(Jdiff,2,1);
    Ji = Jii;
end
figure
loglog(diff_test(1:end-1), variance, ".-", 'LineWidth',1)
legend("$\dot{\lambda_1}$", "$\dot{\lambda_2}$", "$\dot{\lambda_3}$", "$\lambda_1$", "$\lambda_2$", "$\lambda_3$", "$t_f$", 'Interpreter', 'latex', 'Location','southwest')
[~,idx] = min(variance, [], 1, 'omitnan');
% d_lambdaAug = diff_test(idx);
% d_lambdaAug = [0.3, 1.3895e-7, 7.90604e-8, 0.184207, 2.32995e-6, 7.54312e-7, 3.90694e-5];
d_lambdaAug = 1e-6*ones(1,7);
xlabel("Finite Difference Value")
ylabel("Norm of Jacobian Column")
grid on

%% Implement Predictor Corrector

[t_hist, x_hist, orbitEls_hist, A_hist, B_hist, E_hist] = predictorCorrector(x0, maxThrust, maxThrust, [B_guess; A_guess; tf_est], d_lambdaAug, targetOrbitEls, mu);
xf = x_hist(:,end);
% orbitEls_hist(5,:) = orbitEls_hist(4,:) + orbitEls_hist(5,:);
fprintf("Final Position: \n")
disp(xf)
% final_orbit_error = printOrbitError(xf, mu, targetOrbitEls);
final_orbit_error = printOrbitError(orbitEls_hist(:,end), mu, targetOrbitEls);

figure
subplot(2,1,1)
plot(t_hist,A_hist(1,:), "LineWidth",linewidth)
hold on
plot(t_hist,A_hist(2,:), "LineWidth",linewidth)
plot(t_hist,A_hist(3,:), "LineWidth",linewidth)
ylabel("A")
xlabel("time (s)")
legend("$A_1$","$A_2$","$A_3$", 'interpreter', 'latex')
grid on
subplot(2,1,2)
plot(t_hist,B_hist(1,:), "LineWidth",linewidth)
hold on
plot(t_hist,B_hist(2,:), "LineWidth",linewidth)
plot(t_hist,B_hist(3,:), "LineWidth",linewidth)
ylabel("B")
xlabel("time (s)")
legend("$B_1$","$B_2$","$B_3$", 'interpreter', 'latex')
grid on

figure
plot(t_hist, E_hist, 'LineWidth',1)
legend("$\epsilon_1$", "$\epsilon_2$", "$\epsilon_3$", "$\epsilon_4$", "$\epsilon_5$", "$\epsilon_6$", "$\epsilon_7$", 'Interpreter', 'latex', 'Location','southwest')
grid on
xlabel("Time [s]")
ylabel("Error")

figure
orbitLabels = ["$a$","$e$","$i$","$\Omega$","$\omega$","$\nu$"];
desOrbitEls = [a; e; inc; NaN; omega; NaN];
for i=1:6
    subplot(6,1,i)
    plot(t_hist,orbitEls_hist(i,:),"LineWidth",2)
    hold on
    ylabel(orbitLabels(i), 'Interpreter','latex')
    yline(desOrbitEls(i), 'r--')
end
sgtitle("Orbit Elements vs Time")
xlabel("Time [s]")

t_tot = lambdaAug(end);
% A = lambdaAug(4:6);
% B = lambdaAug(1:3);
% [~,x_oneshot] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [0 t_tot], x0, options);
% printOrbitError(x_oneshot(end,:)', mu, targetOrbitEls)

%% Problem 3 Monte Carlo

% runs = 20;
% thrust_nom = maxThrust;
% sigma = 0.05*thrust_nom/3;
% thrust_dispersion = mvnrnd(thrust_nom,sigma^2, runs);  
% % runs = 1;
% % thrust_dispersion = 30.2e-3;
% orbit_error_hist = NaN(4,runs);
% injectionTA = NaN(1,runs);
% burnTimes = NaN(1,runs);
% 
% for i = 1:runs
%     thrust_true = thrust_dispersion(i);
%     lambdaAug_guess = [B_guess; A_guess; tf_est];
%     [t_hist_temp, x_hist_temp, orbitEls_hist_temp, A_hist_temp, B_hist_temp, E_hist_temp] = predictorCorrector(x0, thrust_nom, thrust_true, lambdaAug_guess, d_lambdaAug, targetOrbitEls, mu);
%     [~,orbit_error_hist(:,i)] = printOrbitError(orbitEls_hist_temp(:,end), mu, targetOrbitEls);
%     injectionTA(i) = orbitEls_hist_temp(end,end)*180/pi;
%     burnTimes(i) = t_hist_temp(end);
% end
% 
% data = {orbit_error_hist(1,:), orbit_error_hist(2,:), orbit_error_hist(3,:)*180/pi, orbit_error_hist(4,:)*180/pi, injectionTA, burnTimes};
% labels = {'a Error [km]', 'e Error', 'i Error [degrees]', 'ω Error [degrees]', '𝝂 [degrees]', 'Burn Duration [s]'};
% A = NaN(length(data),2);
% figure
% for i = 1:6
%     [A(i,2),A(i,1)] = std(data{i});
%     subplot(2, 3, i);
%     histogram(data{i}, 15);
%     xlabel(labels{i});
%     ylabel('Frequency');
%     title([labels{i}]);
% end
% sgtitle('Pertubation Histograms');
% 
% meanSTD_table = array2table(A,"VariableNames",["Mean","Std Dev"]);
% meanSTD_table.Properties.RowNames = labels;
% disp(meanSTD_table)

%% Function


function[orbitEls_error, orbitEls_error_primary] = printOrbitError(orbitEls_tf, mu, targetOrbitEls)
% function[orbitEls_error] = printOrbitError(xf, mu, targetOrbitEls)
    % [orbitEls_tf, ~] = rv2orbitEls(xf, mu, 1);
    orbitEls_error = orbitEls_tf - targetOrbitEls;
    combinedOmega = wrapTo2Pi(orbitEls_error(4)+orbitEls_error(5));
    orbitEls_error_primary = [orbitEls_error([1,2,3]); combinedOmega];
    fprintf("Final Orbit Element Errors: \n")
    % disp(orbitEls_error)
    fprintf("a error: %f km\n", orbitEls_error_primary(1))
    fprintf("e error: %f\n", orbitEls_error_primary(2))
    fprintf("i error: %f deg\n", orbitEls_error_primary(3)*180/pi)
    fprintf("omega error: %f deg\n", orbitEls_error_primary(4)*180/pi)
end