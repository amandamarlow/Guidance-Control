function dXdt = lowThrustDynamics(X, L_targ, A_targ, Fmax, k, epsilon, mu)
    r = X(1:3); % extract position from state vector
    v = X(4:6); % extract velocity from state vector
    dr = v; % inertial time derivative of position is simply velocity
    
    F = control_F(r, v, L_targ, A_targ, Fmax, k, epsilon, mu); % get control force
    dv = -mu/(norm(r)^3)*r + F; % inertial acceleration

    dXdt = [dr; dv]; % assign back to state vector derivative
end
