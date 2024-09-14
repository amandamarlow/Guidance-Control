function [nu] = nuFromM(M,a,e,mu)
    % nu = fsolve(@(nu) M-M_from_nu(nu,x), M);
    % [a,e] = unpackState(x);
    E = M2E(M, a, e, mu);
    nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end
