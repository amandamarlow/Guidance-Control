function E = M2E(M, a, e, mu)
    % n = sqrt(mu/a^3);
    E = M;
    
    % Implement Newtons method
    tolerance = 1*10^-6; % rad chosen based on sig figs of information given in question 
    converge = 0;
    for i = 1:1000 % for loop with max iterations to prevent infinite loop if it doesnt converge
        Enext = E - g(E,M,e)/dgdE(E,e);
        if abs(Enext - E) < tolerance
            % fprintf("converged in %d iterations\n",i)
            converge = 1;
            break
        end
        E = Enext;
    end
    if converge == 0
        fprintf("did not converge to tolerance within 1000 iterations")
    end
end

function g = g(E, M, e) 
    % function to find roots of in Newtons method
    g = E - e*sin(E) - M;
end

function dgdE = dgdE(E, e)
    % derivative of root function with respect to E
    dgdE = 1 - e*cos(E);
end