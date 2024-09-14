function [a, e, i, OMEGA, omega, M] = unpackOrbitEls(x)
        a = x(1);
        e = x(2);
        i = x(3);
        omega = x(4);
        OMEGA = x(5);
        M = x(6);
end