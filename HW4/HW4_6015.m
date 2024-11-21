
clear
clc
close all
addpath("C:\Users\marlo\MATLAB Drive\6015")
% constants
mu = 3.986004415e5;


x0 = [8276; 5612; 5; -3.142; 4.672; 0];
% initial guesses
A = [1; 0; 0];
B = [0; 0; 0];
tf = 23; % burnout time
% [~, NO_t0] = rv2orbitEls(x0, mu);
% A_N = NO_t0*A_O;
% B_N = NO_t0*B_O;

a = 10000; % km
i = 0; 
e = 0.0001; % km
omega = 35*pi/180; % deg
targetOrbitEls = [a; e; i; 0; omega; 0];

% [xf, NO_f] = orbitEls2xv(orbitEls, mu);
% [orbitEls0, NO_0] = rv2orbitEls(x0, mu);

rp = a*(1-e);
vp = sqrt(2*mu/rp - mu/a);
rp_O = [rp; 0; 0];
vp_O = [0; vp; 0];
[~, NO_periapsis] = orbitEls2xv([a; e; i; 0; omega; 0], mu);
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
[t,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [0 tf], x0, options);

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
lambdaAug = [B; A; tf];
E = ode45wrapper(x0, lambdaAug, maxThrust, targetOrbitEls, mu);

exponent = linspace(-12, 0, 50);
diff_test = (10*ones(size(exponent))).^exponent;
variance = NaN(size(diff_test, 2)-1, 7);
Ji = jacobian(@(input) ode45wrapper(x0, input, maxThrust, targetOrbitEls, mu), lambdaAug, diff_test(1)*ones(7,1), E);
for i = 1:size(diff_test, 2)-1
    Jii = jacobian(@(input) ode45wrapper(x0, input, maxThrust, targetOrbitEls, mu), lambdaAug, diff_test(i+1)*ones(7,1), E);
    Jdiff = Jii-Ji;
    variance(i,:) = vecnorm(Jdiff,2,1);
    Ji = Jii;
end
figure
loglog(diff_test(1:end-1), variance, ".-", 'LineWidth',1)
legend("$\dot{\lambda_1}$", "$\dot{\lambda_2}$", "$\dot{\lambda_3}$", "$\lambda_1$", "$\lambda_2$", "$\lambda_3$", "$t_f$", 'Interpreter', 'latex', 'Location','southwest')
[~,idx] = min(variance, [], 1, 'omitnan');
% d_lambdaAug = diff_test(idx);
d_lambdaAug = [0.3, 1.3895e-7, 7.90604e-8, 0.184207, 2.32995e-6, 7.54312e-7, 3.90694e-5];

%% Implement Predictor Corrector
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
guidance_updates = 1; % s
% lambdaAug = [B; A; tf];
max_newton_iterations = 10;
max_line_search_iterations = 100;
% error_tol = 1;
error_tol = 0.01;
storeE = [];
Wnorm = @(input) weightedNorm(input, x0, tf);
while tf >= 5
    lambdaAug = [B; A; tf];

    for k = 1:max_newton_iterations
        
        E = ode45wrapper(x0, lambdaAug, maxThrust, targetOrbitEls, mu);
        J = jacobian(@(input) ode45wrapper(x0, input, maxThrust, targetOrbitEls, mu), lambdaAug, d_lambdaAug, E);
        
        gamma = 1;
        for j = 1:max_line_search_iterations+1
            lambdaAug_new = lambdaAug + J\E.*gamma;
            E_new = ode45wrapper(x0, lambdaAug_new, maxThrust, targetOrbitEls, mu);
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

        storeE(:,end+1) = E_new;
        lambdaAug = lambdaAug_new;
        % if norm(E_new) < error_tol
        if failed
            fprintf("line search failed to improve error in %f iterations\n", max_line_search_iterations)
            break
        elseif Wnorm(E_new) < error_tol
            fprintf("Converged! E final: \n")
            disp(norm(E_new))
            disp(Wnorm(E_new))
            break
        elseif k == max_newton_iterations
            fprintf("Failed to converge in %f Newton-Raphson iterations. E final: \n", max_newton_iterations)
            disp(norm(E_new))
            disp(Wnorm(E_new))
        end
        
    end
    tf = lambdaAug(end)-guidance_updates;
    A = lambdaAug(4:6);
    B = lambdaAug(1:3);
    [~,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [0 guidance_updates], x0, options);
    x0 = x(end,:)';
end
[~,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [0 tf], x0, options);
xf = x(end,:)';
fprintf("Final Position: \n")
disp(xf)
[orbitEls_tf, NO_tf] = rv2orbitEls(xf, mu);
orbitEls_error = orbitEls_tf - targetOrbitEls;
orbitEls_error = orbitEls_error([1,2,3,5]);
fprintf("Final Orbit Element Errors: \n")
disp(orbitEls_error)
%%

function[J] = jacobian(fun, lambdaAug, d_lambdaAug, E0)
    J = NaN(7);
    for i = 1:7
        di = zeros(7,1);
        di(i) = d_lambdaAug(i);
        Ei = fun(lambdaAug+di);
        J(:,i) = (E0-Ei)/d_lambdaAug(i);
    end
end

function[Wnorm] = weightedNorm(vec, x0, tf)
    if size(vec,2) >1
        vec = vec';
    end
    W = 10*[1/norm(x0(4:6)).*ones(1,3), 1/norm(x0(1:3)).*ones(1,3), 1/tf];
    W = diag(W);
    Wnorm = sqrt(vec'*W*vec);
end