% function for ode45 to integrate
function dXdt = rdot_rddot(X, mu, perturbing_dv)
    r = X(1:3); % extract position from state vector
    v = X(4:6); % extract velocity from state vector
    dr = v; % inertial time derivative of position is simply velocity
    if exist("perturbing_dv","var")
        dv = -mu/(norm(r)^3)*r + perturbing_dv; % inertial time derivative of velocity from equation 2.22
    else
        dv = -mu/(norm(r)^3)*r; % inertial time derivative of velocity from equation 2.22
    end
    dXdt = [dr; dv]; % assign back to state vector derivative
end

