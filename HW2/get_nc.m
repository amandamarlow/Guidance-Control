function [nc] = get_nc(x, tgo, N, aug)
%NC Summary of this function goes here
%   Detailed explanation goes here

nc = zeros(length(tgo),1);
for i = 1:length(tgo)
    y = x(1,i);
    ydot = x(2,i);
    nt = x(3,i);

    if aug
        ZEM_perp = y + ydot*tgo(i) + nt/2*tgo(i)^2;
    else
        ZEM_perp = y + ydot*tgo(i);
    end
    nc(i) = N*ZEM_perp/tgo(i)^2;
end
end

