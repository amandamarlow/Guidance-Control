function [F, G] = control_F(r, v, L_targ, A_targ, Fmax, k, epsilon, mu)
% control_F takes in the current state and target angular momentum and
% laplace vector along with any tuning parameters and outputs the control
% force vector
%   r = position (must be a matrix of column vectors)
%   v = velocity (must be a matrix of column vectors)
%   L = angular momentum
%   A = Laplace Vector
%   Fmax = maximum thrust level
%   k = tuning parameter
%   epsilon = tuning parameter governing switching condition for control F

% [m,n] = size(r);
% F = NaN(m,n);
% for i = 1:n
    % L = cross(r(:,i), v(:,i));
    % A = cross(v(:,i), L) - mu*r(:,i)/norm(r(:,i));
    % delta_L = L - L_targ;
    % delta_A = A - A_targ;
    % G = -( cross(k*delta_L, r(:,i)) + cross(L, delta_A) + cross(cross(delta_A, v(:,i)), r(:,i)) );
    % if norm(G) < epsilon*Fmax
    %     F(:,i) = G/epsilon;
    % else
    %     F(:,i) = Fmax*G/norm(G);
    % end
% end
    L = cross(r, v);
    A = cross(v, L) - mu*r/norm(r);
    delta_L = L - L_targ;
    delta_A = A - A_targ;
    G = -( cross(k*delta_L, r) + cross(L, delta_A) + cross(cross(delta_A, v), r) );
    if norm(G) < epsilon*Fmax
        F = G/epsilon;
    else
        F = Fmax*G/norm(G);
    end

    if norm(F) < 2e-4
        F = zeros(3,1);
    end
end

