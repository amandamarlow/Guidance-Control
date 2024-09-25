function [t_hist, x_hist] = LambertGuidance(x0_N, rf_N, TOF, at, mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

update_period = 60; % s
dataPointsInPeriod = 10;
tol = 1e-12; % km/s
tspan = 0:update_period:TOF;
% t_hist = zeros(length(tspan)*dataPointsInPeriod,1);
% x_hist = zeros(length(tspan)*dataPointsInPeriod,6);
t_hist = NaN(1);
x_hist = NaN(1,6);
xi = x0_N;
for i = 1:length(tspan)-1
    r_N = xi(1:3);
    % r = norm(r_N);
    v_N = xi(4:6);
    tgo = TOF-tspan(i);
    
    using_long_transfer_angle = 0;
    [vr_N, ~, ~] = Lamberts(mu, r_N, rf_N, using_long_transfer_angle, tgo);
    vg_N = vr_N - v_N;
    vg = norm(vg_N);
    vhat = vg_N/norm(vg_N);
    % if all(vg_N < ones(3,1)*tol)
    % if norm(vg_N) < at*update_period
        % at_N = zeros(3,1);
    % else
        at_N = at*vhat;
    % end
    
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    % [t_hist((i-1)*dataPointsInPeriod+1:i*dataPointsInPeriod),x_hist((i-1)*dataPointsInPeriod+1:i*dataPointsInPeriod,:)] = ode45(@(t,x) rdot_rddot(x, mu, at_N), linspace(tspan(i),tspan(i+1),dataPointsInPeriod), xi, options);
    burnTime = min(vg/at,update_period);
    if burnTime > 1e-6
        [t,x] = ode45(@(t,x) rdot_rddot(x, mu, at_N), [tspan(i) tspan(i)+burnTime], xi, options);
        t_hist(end:end+length(t)-1,1) = t;
        x_hist(end:end+length(t)-1,1:6) = x;
    end
    if burnTime < update_period
        [t,x] = ode45(@(t,x) rdot_rddot(x, mu, zeros(3,1)), [t_hist(end) tspan(i+1)], x(end,:), options);
        t_hist(end:end+length(t)-1,1) = t;
        x_hist(end:end+length(t)-1,1:6) = x;
    end
    xi = x_hist(end,:)';
end

end

