function [dSdt] = xK_ODE(t, S, A, B, Q, R)
%NVM this doesn't work K needs to be integrated backward
[n,m] = size(B);

x = S(1:n);
K = reshape(S(n+1:end),[n,n]);

u = -R\B'*K*x;
dxdt = A*x + B*u;

dKdt = -K*A+K*B/R*B'*K-Q-A'*K;
dKdt = reshape(dKdt,[],1);

dSdt = [dxdt; dKdt];

end

