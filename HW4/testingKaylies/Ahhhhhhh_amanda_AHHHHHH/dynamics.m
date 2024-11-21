function x_dot = dynamics(t,X, mu, B, A)

% acceel
a_fixed = 30E-3; % km/s     % KAYLIE COME BACK PUT INTO CONTROL FUNCTION
u_hat = (A + B.*t)/norm(A + B.*t);
a_vec = a_fixed.*u_hat;

% convert u to inertial by finding true anomoly and w
[~,i,~,w,O,ta] = r_and_v_to_orbital_el(X(1:3)',X(4:6)', mu);
DCM = rth2inertial_DCM(w+ta,i,O);
a_vec = DCM*a_vec; % now in inertial

% update dynamics
x_dot(1:3) = X(4:6);
x_dot(4:6) = a_vec - mu/norm(X(1:3))^3*X(1:3);
x_dot = x_dot';

end