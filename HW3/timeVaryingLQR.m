function [xdot] = xK_ODE(t, x, A, B, R, t_hist, K_hist)
%NVM this doesn't work K needs to be integrated backward

K = squeeze(interp1(t_hist,K_hist,t));
% u = -R\B'*K*x;
u = -inv(R)*B'*K*x;
xdot = A*x + B*u;

end

