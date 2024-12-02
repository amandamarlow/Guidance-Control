function [t_hist, x_hist, orbitEls_hist, A_hist, B_hist, E_hist] = predictorCorrector(x0, thrust_nom, thrust_true, lambdaAug_guess, d_lambdaAug, targetOrbitEls, mu)
%PREDICTORCORRECTOR Summary of this function goes here
%   Detailed explanation goes here

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    guidance_updates = 1; % s
    max_newton_iterations = 10;
    max_line_search_iterations = 20;
    error_tol = 1e-10;
    x_current = x0;
    t_go = lambdaAug_guess(end);
    t = 0;
    A = lambdaAug_guess(4:6);
    B = lambdaAug_guess(1:3);
    A_hist = A;
    B_hist = B;
    t_hist = 0;
    E_hist = [];
    x_hist = x0;
    while t_go >= 2
        lambdaAug = [B; A; t_go];
        tstep = min(guidance_updates, t_go);
    
        for k = 1:max_newton_iterations
            
            E = ode45wrapper(t, x_current, lambdaAug, thrust_nom, targetOrbitEls, mu);
            J = jacobian(@(input) ode45wrapper(t, x_current, input, thrust_nom, targetOrbitEls, mu), lambdaAug, d_lambdaAug, E);
            
            if norm(E) > error_tol
                gamma = 1;
                for j = 1:max_line_search_iterations+1
                    lambdaAug_new = lambdaAug + J\E.*gamma;
                    E_new = ode45wrapper(t, x_current, lambdaAug_new, thrust_nom, targetOrbitEls, mu);
                    if norm(E_new) < norm(E)
                        failed = false;
                        break
                    elseif j == max_line_search_iterations+1
                        E_new = E;
                        lambdaAug_new = lambdaAug;
                        failed = true;
                    end
                    gamma = gamma*0.5;
                end
        
                lambdaAug = lambdaAug_new;
                if failed
                    fprintf("line search failed to improve error on the %d newton-raphson iteration\n", k)
                    disp(norm(E_new))
                    break
                elseif norm(E_new) <= error_tol
                    fprintf("Converged after %d Newton-Raphson iterations! E final: \n", k)
                    disp(norm(E_new))
                    break
                elseif k == max_newton_iterations
                    fprintf("Failed to converge in %f Newton-Raphson iterations. E final: \n", max_newton_iterations)
                    disp(norm(E_new))
                end
            else
                E_new = E;
                % fprintf("Converged before Newton-Raphson \n")
                break
            end
            
        end
        t_go = lambdaAug(end)-tstep;
        A = lambdaAug(4:6);
        B = lambdaAug(1:3);
        [~,x_out] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, thrust_true), [t t+tstep], x_current, options);
        x_current = x_out(end,:)';
        t = t + tstep;
    
        A_hist(:,end+1) = A;
        B_hist(:,end+1) = B;
        t_hist(:,end+1) = t;
        E_hist(:,end+1) = E;
        x_hist(:,end+1) = x_current;
        
    end
    A_hist(:,end+1) = A;
    B_hist(:,end+1) = B;
    t_hist(:,end+1) = t + t_go;
    E_hist(:,end+1) = E_new;
    E_hist(:,end+1) = ode45wrapper(t, x_current, [B;A;t_go], thrust_nom, targetOrbitEls, mu);
    
    [~,x_out] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, thrust_true), [t t+t_go], x_current, options);
    xf = x_out(end,:)';
    x_hist(:,end+1) = xf;
    orbitEls_hist = NaN(6,size(x_hist, 2));
    for i = 1:size(x_hist,2)
        [orbitEls_hist(:,i), ~] = rv2orbitEls(x_hist(:,i), mu, 1);
    end

end

% function[Wnorm] = weightedNorm(vec, x0, tf)
%     if size(vec,2) >1
%         vec = vec';
%     end
%     W = 100*[1/norm(x0(4:6)).*ones(1,3), 1/norm(x0(1:3)).*ones(1,3), 1/tf];
%     W = diag(W);
%     Wnorm = sqrt(vec'*W*vec);
% end