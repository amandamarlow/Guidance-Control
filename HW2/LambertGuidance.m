function [t_hist, x_hist, at_hist] = LambertGuidance(x0_N, rf_N, TOF, at, mu, crossProdSteer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

update_period = 60; % s
tol = 1e-12; % km/s
tspan = 0:update_period:TOF;
rf = norm(rf_N);
G = 3*mu/rf^5*(rf_N*rf_N') - mu/rf^3*eye(3);
t_hist = NaN(1);
x_hist = NaN(1,6);
at_hist = NaN(1);
xi = x0_N;
% for i = 1:length(tspan)-1
t0 = 0;
while t0 < TOF
    tnext = min(t0+update_period, TOF);
    r_N = xi(1:3);
    r = norm(r_N);
    v_N = xi(4:6);
    % tgo = TOF-tspan(i);
    tgo = TOF-t0;
    
    using_long_transfer_angle = 0;
    [vr_N, ~, ~] = Lamberts(mu, r_N, rf_N, using_long_transfer_angle, tgo);
    vg_N = vr_N - v_N;
    vg = norm(vg_N);
    vghat = vg_N/norm(vg_N);
    % if all(vg_N < ones(3,1)*tol)
    % if norm(vg_N) < at*update_period
    %     at_N = zeros(3,1);
    % else
        at_N = at*vghat;
        if exist("crossProdSteer", "var")
            if crossProdSteer == 1
                % G = 3*mu/r^5*(r_N*r_N') - mu/r^3*eye(3);
                Q = -1/(TOF-t0)*(eye(3)+(TOF-t0)^2/2*G);
                b_N = -Q*vg_N;
                b = norm(b_N);
                q = -dot(vghat,b_N) + sqrt(dot(vghat,b_N)^2 - b^2 + at^2);
                at_N = b_N + q*vghat;
            end
        end
    % end

    
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    burnTime = min(vg/at,update_period);
    % burnTime = tnext-t0;
    if burnTime > 1e-9
        % [t,x] = ode45(@(t,x) rdot_rddot(x, mu, at_N), [tspan(i) tspan(i)+burnTime], xi, options);
        [t,x] = ode45(@(t,x) rdot_rddot(x, mu, at_N), [t0 t0+burnTime], xi, options);
        t_hist(end:end+length(t)-1,1) = t;
        x_hist(end:end+length(t)-1,1:6) = x;
        at_hist(end:end+length(t)-1,1) = at;
        % at_hist(end:end+length(t)-1,1) = norm(at_N);
    end
    if burnTime < update_period
        % [t,x] = ode45(@(t,x) rdot_rddot(x, mu, zeros(3,1)), [t_hist(end) tspan(i+1)], x(end,:), options);
        [t,x] = ode45(@(t,x) rdot_rddot(x, mu, zeros(3,1)), [t_hist(end) tnext], x(end,:), options);
        t_hist(end:end+length(t)-1,1) = t;
        x_hist(end:end+length(t)-1,1:6) = x;
        at_hist(end:end+length(t)-1,1) = 0;
    end
    xi = x_hist(end,:)';
    t0 = tnext;
end

end

