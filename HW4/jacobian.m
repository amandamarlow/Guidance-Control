function[J] = jacobian(fun, lambdaAug, d_lambdaAug, E0)
    J = NaN(7);
    for i = 1:7
        di = zeros(7,1);
        di(i) = d_lambdaAug(i);
        Ei = fun(lambdaAug+di);
        J(:,i) = (E0-Ei)/d_lambdaAug(i);
        % J(:,i) = (Ei-E0)/d_lambdaAug(i);
    end
end

