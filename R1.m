function [C] = R1(theta)
%M3 returns the DCM for a rotation of theta about the first axis
% theta must be given in radians

C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
end
