function dXdt = bilinearTangent_dynamics(t, X, mu, A, B, maxThrust)
    r = X(1:3); % extract position from state vector
    v = X(4:6); % extract velocity from state vector
    dr = v; % inertial time derivative of position is simply velocity
    
    u_hat_lvlh = (A+B*t)/norm(A+B*t);
    [~, NO] = rv2orbitEls(X, mu, 1);
    % u_hat_N = R3(pi)*R2(pi/2)*NO*u_hat_O;
    % u_hat_N = NO*R2(-pi/2)*R3(-pi/2)*u_hat_O;
    u_hat_N = NO*R3(-pi/2)*R1(-pi/2)*u_hat_lvlh;
    perturbing_dv = u_hat_N*maxThrust;
    dv = -mu/(norm(r)^3)*r + perturbing_dv; % inertial time derivative of velocity from equation 2.22

    dXdt = [dr; dv]; % assign back to state vector derivative
end

