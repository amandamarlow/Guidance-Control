
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

% [xf, NO_f] = orbitEls2xv(orbitEls, mu);
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



%% 

% lambda_v = A; % assuming t = 0
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% [t,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A_N, B_N, maxThrust), [0 tf], x0, options);
[t,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A_guess, B_guess, maxThrust), [0 tf_est], x0, options);

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

%% Implement Predictor Corrector

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
guidance_updates = 1; % s
% guidance_updates = 0.5; % s
max_newton_iterations = 10;
max_line_search_iterations = 20;
% error_tol = 1;
% error_tol = 0.5;
% error_tol = 0.1;
% error_tol = 0.001;
error_tol = 1e-10;
storeE = [];
x_current = x0;
t_go = tf_est;
Wnorm = @(input) weightedNorm(input, x0, t_go);
t = 0;
A = A_guess;
B = B_guess;
A_hist = A;
B_hist = B;
t_hist = 0;
E_hist = [];
x_hist = x0;
% lambdaAug = [B; A; tf_est];
% for oneIter = 1
while t_go >= 2
    lambdaAug = [B; A; t_go];
    tstep = min(guidance_updates, t_go);

    for k = 1:max_newton_iterations
        
        E = ode45wrapper(t, x_current, lambdaAug, maxThrust, targetOrbitEls, mu);
        J = jacobian(@(input) ode45wrapper(t, x_current, input, maxThrust, targetOrbitEls, mu), lambdaAug, d_lambdaAug, E);
        
        gamma = 1;
        for j = 1:max_line_search_iterations+1
            % lambdaAug_new = lambdaAug - J\E.*gamma;
            lambdaAug_new = lambdaAug + J\E.*gamma;
            E_new = ode45wrapper(t, x_current, lambdaAug_new, maxThrust, targetOrbitEls, mu);
            if norm(E_new) < norm(E)
            % if Wnorm(E_new) < Wnorm(E)
                failed = false;
                break
            elseif j == max_line_search_iterations+1
                % fprintf("failed to improve error in %f gamma iterations\n", max_line_search_iterations)
                E_new = E;
                lambdaAug_new = lambdaAug;
                failed = true;
            end
            gamma = gamma*0.5;
        end

        % storeE(:,end+1) = E_new;
        lambdaAug = lambdaAug_new;
        if failed
            % fprintf("line search failed to improve error in %f iterations\n", max_line_search_iterations)
            fprintf("line search failed to improve error on the %d newton-raphson iteration\n", k)
            disp(norm(E_new))
            disp(Wnorm(E_new))
            break
        elseif norm(E_new) <= error_tol
        % elseif Wnorm(E_new) < error_tol
            fprintf("Converged after %d Newton-Raphson iterations! E final: \n", k)
            disp(norm(E_new))
            disp(Wnorm(E_new))
            break
        elseif k == max_newton_iterations
            fprintf("Failed to converge in %f Newton-Raphson iterations. E final: \n", max_newton_iterations)
            disp(norm(E_new))
            disp(Wnorm(E_new))
        end
        
    end
    t_go = lambdaAug(end)-tstep;
    A = lambdaAug(4:6);
    B = lambdaAug(1:3);
    [~,x_out] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [t t+tstep], x_current, options);
    x_current = x_out(end,:)';
    t = t + tstep;

    A_hist(:,end+1) = A;
    B_hist(:,end+1) = B;
    t_hist(:,end+1) = t;
    E_hist(:,end+1) = E;
    x_hist(:,end+1) = x_current;
    
end
A_hist(:,end+1) = A;
B_hist(:,end+1) = B;
t_hist(:,end+1) = t + t_go;
E_hist(:,end+1) = E_new;
% lambdaAug(end) = tf;
E_hist(:,end+1) = ode45wrapper(t, x_current, [B;A;t_go], maxThrust, targetOrbitEls, mu);

[~,x_out] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [t t+t_go], x_current, options);
xf = x_out(end,:)';
x_hist(:,end+1) = xf;
orbitEls_hist = NaN(6,size(x_hist, 2));
for i = 1:size(x_hist,2)
    [orbitEls_hist(:,i), ~] = rv2orbitEls(x_hist(:,i), mu, 1);
end
fprintf("Final Position: \n")
disp(xf)
printOrbitError(xf, mu, targetOrbitEls)

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

figure
ylabels = ["a","e","i","W","w","ta"];
desOrbitEls = [a; e; inc; NaN; omega; NaN];
for i=1:6
    subplot(6,1,i)
    plot(t_hist,orbitEls_hist(i,:),"LineWidth",2)
    hold on
    ylabel(ylabels(i))
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

runs = 20;
thrust_nom = maxThrust;
sigma = 0.05*thrust_nom/3;
thrust_dispersion = mvnrnd(thrust_nom,sigma^2, runs);              

%% Function

function[J] = jacobian(fun, lambdaAug, d_lambdaAug, E0)
    J = NaN(7);
    for i = 1:7
        di = zeros(7,1);
        di(i) = d_lambdaAug(i);
        Ei = fun(lambdaAug+di);
        J(:,i) = (E0-Ei)/d_lambdaAug(i);
        % J(:,i) = (Ei-E0)/d_lambdaAug(i);
    end
end

function[Wnorm] = weightedNorm(vec, x0, tf)
    if size(vec,2) >1
        vec = vec';
    end
    W = 100*[1/norm(x0(4:6)).*ones(1,3), 1/norm(x0(1:3)).*ones(1,3), 1/tf];
    W = diag(W);
    Wnorm = sqrt(vec'*W*vec);
end

function[] = printOrbitError(xf, mu, targetOrbitEls)
    [orbitEls_tf, ~] = rv2orbitEls(xf, mu, 1);
    orbitEls_error = orbitEls_tf - targetOrbitEls;
    orbitEls_error = orbitEls_error([1,2,3,5]);
    fprintf("Final Orbit Element Errors: \n")
    % disp(orbitEls_error)
    fprintf("a error: %f km\n", orbitEls_error(1))
    fprintf("e error: %f\n", orbitEls_error(2))
    fprintf("i error: %f deg\n", orbitEls_error(3)*180/pi)
    fprintf("omega error: %f deg\n", orbitEls_error(4)*180/pi)
end