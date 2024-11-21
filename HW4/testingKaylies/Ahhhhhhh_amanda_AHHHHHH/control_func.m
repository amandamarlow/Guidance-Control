function a_vec_matrix = control_func(t,X, mu, A, B)

% % acceel
% a_fixed = 30E-3; % km/s     % KAYLIE COME BACK PUT INTO CONTROL FUNCTION
% u_hat = (A + B.*t)/norm(A + B.*t);
% a_vec = a_fixed.*u_hat;
% 
% % convert u to inertial by finding true anomoly and w
% [~,i,~,w,O,ta] = r_and_v_to_orbital_el(X(1:3)',X(4:6)', mu);
% DCM = rth2inertial_DCM(w+ta,i,O);
% a_vec = DCM*a_vec; % now in inertial

 % Ensure t is a column vector
    t = t(:);

    % Fixed acceleration
    a_fixed = 30E-3; % km/s

    % Calculate unit vectors u_hat for each time step
    u = A + B .* t; % This will automatically expand t to match A and B dimensions
    u_norm = sqrt(sum(u.^2, 2)); % Compute the norm for each row
    u_hat = u ./ u_norm; % Normalize each row

    % Calculate acceleration vectors
    a_vec = a_fixed * u_hat;

    % Convert to inertial frame
    [~, i, ~, w, O, ta] = r_and_v_to_orbital_el(X(:,1:3), X(:,4:6), mu);
    DCM = rth2inertial_DCM(w + ta, i, O);

    % Apply DCM to all acceleration vectors
    a_matrix = (DCM * a_vec')'; % Transpose to apply DCM and then transpose 

end