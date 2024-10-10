function [dKdt] = K_ODE(t, K, A, B, Q, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

K = reshape(K,5,5);

dKdt = -K*A+K*B/R*B'*K-Q-A'*K;
dKdt = reshpape(dKdt,[],1);
end

