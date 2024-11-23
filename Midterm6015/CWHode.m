function [dXdt] = CWHode(t, X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 3.986004418e14; % m^3/s^2
a = 6728e3; % m
n = sqrt(mu/a^3);

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    3*n^2, 0, 0, 2*n;
    0, 0, -2*n, 0];
B = [0, 0;
    0, 0;
    1, 0;
    0, 1];

u = zeros(2,1);

dXdt = A*X + B*u;

end

