function weighted_norm = weighted_norm_func(E_vec)
% function to more cleanly modify the weight of each parameter, these
% weights were determined in the expiement script

% weighting
% e1w = 1813599.25900668*E_vec(1);
% e2w = 397.069173427775*E_vec(2);
% e3w = 221.823568429933*E_vec(3);
% e4w = 0.111135840858210*E_vec(4);
% e5w = 0.199128609795830*E_vec(5);
% e6w = 96.3135671471813*E_vec(6);
% e7w = 42.4396068662541*E_vec(7);

e1w = E_vec(1);
e2w = E_vec(2);
e3w = E_vec(3);
e4w = E_vec(4);
e5w = E_vec(5);
e6w = E_vec(6);
e7w = E_vec(7);

E_weighted = [e1w e2w e3w e4w e5w e6w e7w];

weighted_norm = norm(E_weighted);

% function[Wnorm] = weightedNorm(vec, x0, tf)
%     if size(vec,2) >1
%         vec = vec';
%     end
%     W = 10*[1/norm(x0(4:6)).*ones(1,3), 1/norm(x0(1:3)).*ones(1,3), 1/tf];
%     W = diag(W);
%     Wnorm = sqrt(vec'*W*vec);
% end

end