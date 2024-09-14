% function for ode45 to integrate
function dXdt = rdot_rddot(X, mu)
    r = X(1:3); % extract position from state vector
    v = X(4:6); % extract velocity from state vector
    dr = v; % inertial time derivative of position is simply velocity
    dv = -mu/(norm(r)^3)*r; % inertial time derivative of velocity from equation 2.22
    dXdt = [dr; dv]; % assign back to state vector derivative
end

