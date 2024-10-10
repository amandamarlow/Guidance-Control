function [dpdt] = adjointODE(t, p, A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dpdt = -A'*p;
end

