function [r, v] = orbital_elements_to_state_vectors(a, i, e, w, RAAN, nu)
% orbital_elements_to_state_vectors - Convert orbital elements to ECI state vectors
%
% Inputs:
%   a    - Semi-major axis (km)
%   i    - Inclination (degrees)
%   e    - Eccentricity
%   w    - Argument of perigee (degrees)
%   RAAN - Right ascension of ascending node (degrees)
%   nu   - True anomaly (degrees)
%
% Outputs:
%   r    - Position vector in ECI frame (km)
%   v    - Velocity vector in ECI frame (km/s)

% Earth's gravitational parameter (km^3/s^2)
mu = 398600.4418;

% Convert angles to radians
i = deg2rad(i);
w = deg2rad(w);
RAAN = deg2rad(RAAN);
nu = deg2rad(nu);

% Calculate auxiliary parameters
p = a * (1 - e^2);
h = sqrt(mu * p);

% Compute position and velocity in perifocal frame
r_w = (h^2 / mu) * (1 / (1 + e * cos(nu))) * [cos(nu); sin(nu); 0];
v_w = (mu / h) * [-sin(nu); e + cos(nu); 0];

% Compute rotation matrix from perifocal to ECI frame
R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
R3_w = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];
Q_wx = R3_W * R1_i * R3_w;

% Transform position and velocity to ECI frame
r = Q_wx * r_w;
v = Q_wx * v_w;

% Ensure column vectors
r = r(:);
v = v(:);

end