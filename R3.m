function [C] = R3(theta)
%M3 returns the DCM for a rotation of theta about the third axis
% theta must be given in radians

C = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
end



