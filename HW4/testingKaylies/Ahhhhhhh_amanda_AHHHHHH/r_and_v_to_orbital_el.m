function [a,i,e,w,Omega,theta_star] = r_and_v_to_orbital_el(r_vec,v_vec, mu)
% converting from r and v vetcors (in inertial) to e,a,i, and w
% r and v are n by 3

% Finding h
h_vec = cross(r_vec, v_vec, 2);
h_mag = vecnorm(h_vec, 2, 2);

% Finding e
e_vec = cross(v_vec, h_vec, 2) ./ mu - r_vec ./ vecnorm(r_vec, 2, 2);
e = vecnorm(e_vec, 2, 2);

% Finding a
a = h_mag.^2 ./ (mu * (1 - e.^2));

% Finding i
z_hat = [zeros(size(r_vec,1),1) zeros(size(r_vec,1),1) ones(size(r_vec,1),1)];
i = acosd(dot(h_vec, z_hat, 2) ./ h_mag);

% Finding w (little omega)
n_vec = cross(z_hat, h_vec, 2);
n_mag = vecnorm(n_vec, 2, 2);
w = acosd(dot(n_vec, e_vec, 2) ./ (e .* n_mag));
% Adjust w when e_vec.z < 0
w(dot(e_vec, z_hat, 2) < 0) = -w(dot(e_vec, z_hat, 2) < 0); % check that this is valid

% finding RAAN ;
x_hat = [ones(size(r_vec,1),1) zeros(size(r_vec,1),1) zeros(size(r_vec,1),1)];
y_hat = [zeros(size(r_vec,1),1) ones(size(r_vec,1),1) zeros(size(r_vec,1),1)];
Omega = acosd(dot(n_vec, x_hat,2)./n_mag);
% quad check
Omega(dot(n_vec, y_hat, 2) < 0) = -Omega(dot(n_vec, y_hat, 2) < 0); % check that this is valid

% true anamoly
theta_star = acosd(dot(r_vec,e_vec, 2)./(vecnorm(r_vec, 2, 2).*e));
theta_star(dot(r_vec, v_vec, 2) < 0) = -theta_star(dot(r_vec, v_vec, 2) < 0); % check that this is valid
% if dot(r1,v1) < 0
%     theta_star = -theta_star; % assuming degrees
% end


end