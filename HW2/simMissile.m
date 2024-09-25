function [t,x,nc] = simMissile(R0, Vm, Vt, N, HE0, nt, aug)
%SIMMISSILE Summary of this function goes here
%   Detailed explanation goes here

ydot0 = Vm*tand(HE0);
x0 = [0; ydot0; nt];

Vc = Vm + Vt;
tgo = R0/Vc;
tspan = [0, tgo];
[t,x] = ode45(@(t,x) proNavLinear(t, x, N, tgo, aug), tspan, x0);
tgo_vec = tgo-t';
nc = get_nc(x', tgo_vec, N, aug);


if aug
    leg = strcat("Aug N' = ", num2str(N));
    lineStyle = "--";
else
    leg = strcat("N' = ", num2str(N));
    lineStyle = "-";
end

% labels = ["$y$ (ft)", "$\dot{y}$ (ft/s)", "$n_t$ (ft/s$^2$)"];
% for i = 1:3
%     subplot(3,1,i)
%     hold on
%     plot(t, x(:,i),lineStyle, "DisplayName", leg)
%     ylabel(labels(i), 'Interpreter', 'latex')
%     legend('Location', 'best')
% end

labels = ["$y$ (ft)", "$n_c$ (ft/s$^2$)"];
subplot(2,1,1)
hold on
plot(t, x(:,1), lineStyle, "DisplayName", leg)
ylabel(labels(1), 'Interpreter', 'latex')
legend('Location', 'best')
subplot(2,1,2)
hold on
plot(t, nc, lineStyle, "DisplayName", leg)
ylabel(labels(2), 'Interpreter', 'latex')
legend('Location', 'best')

xlabel("time (s)")
end

