function DCM_i_to_rot = R3_DCM(angle)
% angle is w PLUS true anamoly (in degrees)
% this is the rotation from inertial to rot frame

DCM_i_to_rot = [cosd(angle) sind(angle) 0 ; -sind(angle) cosd(angle) 0 ; 0 0 1];


end