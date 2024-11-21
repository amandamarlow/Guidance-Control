function DCM = rth2inertial_DCM(o,i,t)
% convert to rads
o = deg2rad(o); % big omega
i = deg2rad(i); % inclination
t = deg2rad(t); % ta plus lil omega

DCM = [cos(o)*cos(t)-sin(o)*cos(i)*sin(t) -cos(o)*sin(t)-sin(o)*cos(i)*cos(t) sin(o)*sin(i);...
    sin(o)*cos(t)+cos(o)*cos(i)*sin(t) -sin(o)*sin(t)+cos(o)*cos(i)*cos(t) -cos(o)*sin(i);...
    sin(i)*sin(t) sin(i)*cos(t) cos(i)];

end