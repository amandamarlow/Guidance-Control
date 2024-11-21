%% experimenting around with weightednorm stuff
% desired final state
des_orbital_el = [10000 0 0.0001 35 0]; % a [km], i, e, w [deg], big O [deg]
mu = 3.986E5; % earthtf = 23;
X_0 = [8276; 5612; 5; -3.142; 4.672; 0]; % km and km/s
tf = 23;
% inital control
A = [1 0 0]';
B = [0 0 0]';
u_vec = [B; A; tf]; % initalizing control vector
t_span = [0 u_vec(7)];             
[t_nom, X_nom] = ode45(@(t_nom,X_nom) dynamics(t_nom ,X_nom ,mu, u_vec(1:3), u_vec(4:6)), t_span, X_0);
E0 = E_func(X_nom(end,:), tf, mu, des_orbital_el, u_vec);

% [r, v] = orbital_elements_to_state_vectors(10001, 0.01, 0.00011, 35.01, 0.01, 0.01);
% [r, v] = orbital_elements_to_state_vectors(10010, 0.05, 0.0002, 35.05, 0.05, 0.05);
% [r, v] = orbital_elements_to_state_vectors(10000, 0.01, 0.0001, 35, 0, 0);
[r, v] = orbital_elements_to_state_vectors(10020, 0.1, 0.00025, 35.1, 0.1, 0.1);
state_tf = [r;v]';
A = [1 0.001 0.001]';
B = [0.001 0.001 0.001]';
u_vec = [B; A; tf]; % initalizing control vector
E = E_func(state_tf, tf, mu, des_orbital_el, u_vec)
weighted_norm = weighted_norm_func(E)
%  1.0e+06 *
% 
%    1.813599259006683
%    0.000397069173428
%   -0.000221823568430
%   -0.000000111135841
%    0.000000199128610
%   -0.000251279538342
%                  Inf
% 
%   1.0e+06 *
% 
%    1.813599259006683
%    0.000397069173428
%   -0.000221823568430
%   -0.000000111135841
%    0.000000199128610
%   -0.000096313567147
%    0.000042439606866



%% J convergence study
des_orbital_el = [10000 0 0.0001 35 0]; % a [km], i, e, w [deg], big O [deg]
mu = 3.986E5; % earthtf = 23;
X_0 = [8276; 5612; 5; -3.142; 4.672; 0]; % km and km/s
tf = 23;
time_vec = [];
time_vec = [time_vec 0]; % kaylie come back and implement
% inital control
A = [1 0 0]';
B = [0 0 0]';
% convert to inertial
% get orbital elements
% [~,i,~,w,O,ta] = r_and_v_to_orbital_el(X_0(1:3)',X_0(4:6)', mu);
% DCM_i_to_rot = rth2inertial_DCM(w+ta,i,O);
% A = DCM_i_to_rot'*A; % now in inertial
% B = DCM_i_to_rot'*B; % now in inertial

threshold = 1; % Kaylie come back and change
u_vec = [B; A; tf]; % initalizing control vector

%E0 = [-14322.0384894209 3.12812153156373 -5.46839848641512 -4332.31860268216 8272.40093474971 -0.00194072532411303 0];
E0 = [-0.692822992207353 0.688519859266402 0.0000529806032  -8.553230318962960  5.003440557648769 -0.003979985378460      0];
exponent = linspace(-15, 0, 30);
%exponent = linspace(-5, 0, 30);
diff_vec = (10*ones(size(exponent))).^exponent;
J_norm = zeros(7,length(diff_vec));
for i = 1:length(diff_vec)
    J_norm(:,i) = J_func_study(E0', u_vec, mu, des_orbital_el, X_0, diff_vec(i));
end

J_diff = abs(diff(J_norm,1,2));

figure
for i= 1:7
    loglog(diff_vec(1:end-1),J_diff(i,:))
    hold on
end
legend('1','2','3','4','5','6','7')
xlabel('pertubation')
ylabel('Delta J')

figure
for i= 1:7
    loglog(diff_vec(2:end),J_diff(i,:))
    hold on
end

function J_norm = J_func_study(E0, u_vec, mu, des_orbital_el, X_0, diff_vec)
% getting the Jacobian based on current inital state and final state and
% control

% Taking things out of U_vec
A = u_vec(1:3);
B = u_vec(4:6);
tf = u_vec(7);
t_span = [0 tf];

J = zeros(7,7);

% A
for i = 1:3
    A_i_og = A(i);
    A(i) = A(i) + diff_vec;
    [t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, A, B), t_span, X_0);
    u_vec = [A; B; tf]; % change later but have to send it in right
    E_i = E_func(X_J(end,:), tf, mu, des_orbital_el,u_vec);
    J(:,i) = (E_i-E0)/diff_vec;
    A(i) = A_i_og; % reverting it back
end

% B
for i = 1:3
    B_i_og = B(i);
    B(i) = B(i) + diff_vec;
    [t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, A, B), t_span, X_0);
    u_vec = [A; B; tf]; % change later but have to send it in right
    E_i = E_func(X_J(end,:), tf, mu, des_orbital_el,u_vec);
    J(:,i+3) = (E_i-E0)/diff_vec;
    B(i) = B_i_og; % reverting it back
end

% tf
tf_pert = tf + diff_vec;
t_span = [0 tf_pert];
[t, X_J] = ode45(@(t,X_J) dynamics(t ,X_J ,mu, A, B), t_span, X_0);
u_vec = [A; B; tf_pert];
E_i = E_func(X_J(end,:), tf_pert, mu, des_orbital_el,u_vec);
J(:,7) = (E_i-E0)/diff_vec;

J_norm = vecnorm(J,2,1);

end