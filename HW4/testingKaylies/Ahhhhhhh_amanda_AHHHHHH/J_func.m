function J = J_func(E0, u_vec, mu, des_orbital_el, X_0)
% getting the Jacobian based on current inital state and final state and
% control

% Taking things out of U_vec
B = u_vec(1:3);
A = u_vec(4:6);
tf = u_vec(7);
t_span = [0 tf];

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%diff_vec = 1E-4; % kaylie come back and change likley needs to be different for each one
%diff_vec = [0.3 1.9E-7 2.2E-5 0.1 5.7E-8 2.1E-6 2.7E-7];
diff_vec = [0.3, 1.3895e-7, 7.90604e-8, 0.184207, 2.32995e-6, 7.54312e-7, 3.90694e-5];
% diff_vec = [1e-1, 1e-9,1e-9,1e-1,1e-6,1e-6,1e-5];
J = zeros(7,7);

% A
for i = 1:3
    A_i_og = A(i);
    A(i) = A(i) + diff_vec(i+3);
    [t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, B, A), t_span, X_0, options);
    u_vec = [B; A; tf]; % change later but have to send it in right
    E_i = E_func(X_J(end,:), tf, mu, des_orbital_el,u_vec);
    J(:,i+3) = (E_i-E0)/ diff_vec(i+3);
    A(i) = A_i_og; % reverting it back
end

% B
for i = 1:3
    B_i_og = B(i);
    B(i) = B(i) + diff_vec(i);
    [t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, B, A), t_span, X_0, options);
    u_vec = [B; A; tf]; % change later but have to send it in right
    E_i = E_func(X_J(end,:), tf, mu, des_orbital_el,u_vec);
    J(:,i) = (E_i-E0)/diff_vec(i);
    B(i) = B_i_og; % reverting it back
end

% tf
tf_pert = tf + diff_vec(7);
t_span = [0 tf_pert];
[t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, B, A), t_span, X_0, options);
u_vec = [B; A; tf_pert];
E_i = E_func(X_J(end,:), tf_pert, mu, des_orbital_el,u_vec);
J(:,7) = (E_i-E0)/diff_vec(7);

end