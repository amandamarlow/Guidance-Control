% Kaylie Rick 6015 Hw4 code
clc
clear
close all

%% Problem 1
A = [1 0 0]';
B = [0 0 0]';

X_0 = [8276; 5612; 5; -3.142; 4.672; 0]; % km and km/s
t_span = [0 23]; % seconds

% ode to prop forward
mu = 3.986004415E5; % km units
[t, X] = ode45(@(t,X) dynamics(t,X,mu, B, A), t_span, X_0);

% getting orbital elements
[a,i,e,w,ta] = r_and_v_to_orbital_el(X(:,1:3),X(:,4:6), mu);
% 
% % plots
% figure
% title('Spacecraft Trajectory')
% plot3(X(:,1),X(:,2),X(:,3))
% xlabel('X [km]')
% ylabel('Y [km]')
% zlabel('Z [km]')
% 
% figure
% subplot(1,2,1)
% plot(t,X(:,1:3))
% legend('X','Y','Z')
% subplot(1,2,2)
% plot(t,X(:,4:6))
% legend('VX','VY','VZ')
% 
% figure
% subplot(2,2,1)
% plot(t,a)
% title('a (semi-major axis) over Time')
% xlabel('Time [s]')
% ylabel('a [km]')
% subplot(2,2,2)
% plot(t,e)
% title('e (Eccentricity) over Time')
% xlabel('Time [s]')
% ylabel('e [unitless]')
% subplot(2,2,3)
% plot(t,i)
% title('i (Inclination) over Time')
% xlabel('Time [s]')
% ylabel('i [Degrees]')
% subplot(2,2,4)
% plot(t,w)
% title('w (Arugment of Periapsis) over Time')
% xlabel('Time [s]')
% ylabel('w [Degrees]')

%% Problem 2
% desired final state
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

% big state over time matrix
all_states_big_matrix_yay = [];
all_control_big_matrix_yay = [];
all_time = [];
all_x_finals = [];

delta_t = 0.05; % KAYLIE COME BACK AND CHANGE
threshold = 5; % Kaylie come back and change
u_vec = [B; A; tf]; % initalizing control vector
iterations = 0;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
 weighted_norm_E0 = 100;
while (u_vec(7) > 10) % change to 0 or 5 or whatever.
 
    iterations_comp_100 = [0 1E10];
    while (weighted_norm_E0 > threshold && (abs(iterations_comp_100(1)-iterations_comp_100(2)) > 1)) % add in dynamics for orr it stops changing kaylie come back
        fprintf('Weighted norm = %f and it = %f and time final = %f \n',weighted_norm_E0, iterations, u_vec(7))
        if (iterations == 1 || mod(iterations,10) == 0) % KAYLIE COME BACK AND MOD WITH 100
            iterations_comp_100 = [iterations_comp_100(2) weighted_norm_E0];
        end
        iterations = iterations + 1;
    
        % integrate nominal to current final time
        t_span = [0 u_vec(7)];             
        [t_nom, X_nom] = ode45(@(t_nom,X_nom) dynamics(t_nom ,X_nom ,mu, u_vec(1:3), u_vec(4:6)), t_span, X_0, options);
    
        % compute Jacobian (forward finite differencing)
        % getting E0
        % Getting E0 
        E0 = E_func(X_nom(end,:), u_vec(7), mu, des_orbital_el, u_vec);
%         J = J_func(E0, u_vec, mu, des_orbital_el, X_0); % got rid of outer time
        J = jacobian(@(input) dynamics(x0, input, 30E-3, targetOrbitEls, mu), lambdaAug, diff_test(1)*ones(7,1), E);

        % Line search (inner loop)
        weighted_E0_new_norm = 100E8; % random high number for first check
        weighted_E0_norm_ls = weighted_norm_func(E0); % KAYLIE COME BACK AND MAKE A WEIGHTED NORM FUNC
        gamma = 1;
        while weighted_E0_new_norm > weighted_E0_norm_ls
            u_vec_new = u_vec - J\(gamma*E0);
            % update dynamics with new u
            [t_ls, X_ls] = ode45(@(t_ls,X_ls) dynamics(t_ls ,X_ls ,mu, u_vec_new(1:3), u_vec_new(4:6)), [0 u_vec_new(7)], X_0, options);
            E0_new = E_func(X_ls(end,:), u_vec_new(7), mu, des_orbital_el, u_vec_new);
            weighted_E0_new_norm = weighted_norm_func(E0_new);
            gamma = 0.5*gamma; % kaylie come back and mess around with 
        end % if while loop breaks than this is new E0 for middle loop (update u_vec)
        %fprintf("diff in E_0 = %f \n",E0_new - E0)
        fprintf('Weighted norm = %f and it = %f and time final = %f \n',weighted_norm_E0, iterations, u_vec(7))
        E0 = E0_new;
        % weighted norm E0
        weighted_norm_E0 = weighted_norm_func(E0);
        u_vec = u_vec_new;
    end % once outter loop is satified, update state by running again
    
    % taking an actual time step forward KAYLIE COME BACK AND SAVE
    [t_up, X_up] = ode45(@(t_up,X_up) dynamics(t_up ,X_up ,mu, u_vec_new(1:3), u_vec_new(4:6)), [0 delta_t], X_0,options);
    X_0 = X_up(end,:); % idk if these are right dimentions KAYLIE COME BACK
    u_vec(7) = u_vec(7) - delta_t;

    % save stuff
    all_x_finals = [all_x_finals; X_up(end,:)];
%     all_time = [all_time; time_to_add];
%     all_states_big_matrix_yay = [all_states_big_matrix_yay; X_up]; % KAYLIE COME BACK AND MAKE DIFFERENT TO SAVE ON SPACE
%     all_control_big_matrix_yay = [all_control_big_matrix_yay; u_vec_new'];
end

%% post calcs
% Getting orbial elements at all timesteps
%[a,i,e,w,ta] = r_and_v_to_orbital_el(X(:,1:3),X(:,4:6), mu);
[a_vec,i_vec,e_vec,w_vec,Omega_vec,theta_star_vec] = r_and_v_to_orbital_el(all_x_finals(:,1:3),all_x_finals(:,4:6), mu);
% geting time
data_length = length(all_x_finals);
time_vec = 0:delta_t:(data_length-1)*delta_t;
time_indicies = find(ismember(all_time, time_vec));

% getting control
% a_vec_matrix = control_func(time_vec,all_states_big_matrix_yay, mu, all_control_big_matrix_yay(:,1:3), all_control_big_matrix_yay(:,4:6));

% Plots
% states
figure
subplot(2,1,1)
sgtitle('Inertial State Over Time')
plot(time_vec,all_x_finals(:,1:3))
legend('X Pos','Y Pos','Z Pos')
xlabel('Time [s]')
ylabel('Position [km]')
subplot(2,1,2)
plot(time_vec,all_x_finals(:,4:6))
legend('X Vel','Y Vel','Z Vel')
xlabel('Time [s]')
ylabel('Velocity [km]')

% orbital elements
des_orbital_el = [des_orbital_el NaN];
orbital_elements = {a_vec, i_vec, e_vec, w_vec, Omega_vec, theta_star_vec};
element_labels = {'Semi-major Axis [units]', 'Inclination [degrees]', 'Eccentricity', ...
                  'Argument of Periapsis [degrees]', 'Longitude of Ascending Node [degrees]', ...
                  'True Anomaly [degrees]'};
figure
for i = 1:length(orbital_elements)
    subplot(length(orbital_elements), 1, i);  % Create subplot for each element
    plot(time_vec, orbital_elements{i}, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(des_orbital_el(i), 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel(element_labels{i});  % Customize this as needed
%     title(['Orbital Element: ', element_labels{i}]);
    hold off;
end
sgtitle('Orbital Elements Over Time');


% Plots
% states
figure
subplot(2,1,1)
sgtitle('Inertial State Over Time')
plot(all_time,all_states_big_matrix_yay(:,1:3))
legend('X Pos','Y Pos','Z Pos')
xlabel('Time [s]')
ylabel('Position [km]')
subplot(2,1,2)
plot(all_time,all_states_big_matrix_yay(:,4:6))
legend('X Vel','Y Vel','Z Vel')
xlabel('Time [s]')
ylabel('Velocity [km]')

% control plot
figure
subplot(3,1,1)
sgtitle('LVLH Control Over Time')
plot(time_vec,all_control_big_matrix_yay(:,1:3))
legend('A1','A2','A3')
xlabel('Time [s]')
subplot(3,1,2)
plot(time_vec,all_control_big_matrix_yay(:,4:6))
legend('B1','B2','B3')
xlabel('Time [s]')
subplot(3,1,3)
% plot(all_time,a_vec_matrix) % NOT WORKING RIGHT NOW



% orbital elements
des_orbital_el = [des_orbital_el NaN];
orbital_elements = {a_vec, i_vec, e_vec, w_vec, Omega_vec, theta_star_vec};
element_labels = {'Semi-major Axis [units]', 'Inclination [degrees]', 'Eccentricity', ...
                  'Argument of Periapsis [degrees]', 'Longitude of Ascending Node [degrees]', ...
                  'True Anomaly [degrees]'};
figure
for i = 1:length(orbital_elements)
    subplot(length(orbital_elements), 1, i);  % Create subplot for each element
    plot(all_time, orbital_elements{i}, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(des_orbital_el(i), 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel(element_labels{i});  % Customize this as needed
%     title(['Orbital Element: ', element_labels{i}]);
    hold off;
end
sgtitle('Orbital Elements Over Time');
