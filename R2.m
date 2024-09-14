function [C] = R2(theta)
%M2 returns the DCM for a rotation of theta about the second axis
% theta must be given in radians

C = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
end


