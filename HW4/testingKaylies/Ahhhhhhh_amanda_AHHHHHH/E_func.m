function E = E_func(state_tf, tf, mu, des_orbital_el, u_vec)
% everything in here is in the 

B = u_vec(1:3);
A = u_vec(4:6);
tf = u_vec(7);

% getting final desired orbital elements out and vecs
a_des = des_orbital_el(1);
e_des = des_orbital_el(3);
w_des = des_orbital_el(4);
i_des = des_orbital_el(2);
O_des = des_orbital_el(5);

% convert to current LVLH frame
r_vec = state_tf(1:3);
v_vec = state_tf(4:6);
[~,i_tf,~,w_tf,O_tf,ta_tf] = r_and_v_to_orbital_el(r_vec,v_vec, mu);
inertial_to_lvlh = rth2inertial_DCM(ta_tf+w_tf,i_tf,O_tf)';
r_vec_lvlh = inertial_to_lvlh*r_vec';
v_vec_lvlh = inertial_to_lvlh*v_vec';
tf  = u_vec(7);

% getting rp and vp in periapsis (desired) LVLH frame
% getting them in rotated frame
rp = a_des*(1-e_des);
rp_other_lvlh = [rp 0 0]';
% DCM_i_to_rot = R3_DCM(w);
% rp_inertial = DCM_i_to_rot'*rp_rot;
vp = sqrt(2*mu/rp - mu/a_des);
vp_other_lvlh = [0 vp 0]';
% vp_inertial = DCM_i_to_rot'*vp_rot;

% rotate p's to correct frame
other_lvlh_to_inertial = rth2inertial_DCM(w_des,i_des,O_des);
rp_lvlh = inertial_to_lvlh*other_lvlh_to_inertial*rp_other_lvlh;
vp_lvlh = inertial_to_lvlh*other_lvlh_to_inertial*vp_other_lvlh;


% true anamoly
rp_lvlh_unit = rp_lvlh./norm(rp_lvlh);
vp_lvlh_unit = vp_lvlh./norm(vp_lvlh);
ta = atan2(dot(r_vec_lvlh,vp_lvlh_unit),dot(r_vec_lvlh,rp_lvlh_unit));

% p
p = vp^2*rp^2/mu;

% Getting the first three episolons from v's
VR = -1/rp*sqrt(mu/p)*sin(ta).*rp_lvlh + (1-rp/p*(1-cos(ta))).*vp_lvlh;
e1_to_3 = VR' - v_vec_lvlh'; % make sure correct dimentions KAYLIE COME BACK

% getting r for e4
RR = p/(1+e_des*cos(ta)); % think this is des KAYLIE COME BACK
e4 = RR - norm(r_vec_lvlh);

% getting H and e
H_r = cross(rp_lvlh,VR); % KAYLIE COME BACK REALLY IFY ABOUT THIS
e5 = dot(H_r./norm(H_r),r_vec_lvlh);

% lambda vec
lambda_vec = A + B.*tf;

% e6
g_vec = -mu/norm(r_vec_lvlh)^3.*r_vec_lvlh;
e6 = dot(lambda_vec,g_vec) - dot(B,v_vec_lvlh);

% e7
e7 = sqrt(dot(lambda_vec,lambda_vec)) - 1;

% outputing as E
E = [e1_to_3, e4, e5, e6, e7]'; % KAYLIE COME BACK AND SEE IF THIS WORKS


end