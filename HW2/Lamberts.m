function [v1_N, v2_N, iterations] = Lamberts(mu, r1_N, r2_N, using_long_transfer_angle, TOF)
% function [a, p] = Lamberts(mu, R1_SCI, R2_SCI, using_long_transfer_angle, TOF)
    R1 = norm(r1_N);
    R2 = norm(r2_N);
    % space triangle
    [c,s] = getChordAndSemiPerimeter(r1_N, r2_N, using_long_transfer_angle);
    
    % semi major axis of minimum energy ellipse
    a_min = s/2;
    % calculate time of flight for minimum energy ellipse
    n_min = sqrt(mu/a_min^3);
    [alpha_min, beta_min] = getAlphaBeta(c, s, a_min, 0, using_long_transfer_angle);
    TOF_min = 1/n_min*(alpha_min-beta_min-(sin(alpha_min)-sin(beta_min)));
    % booleon to show if we are using the longer time of flight
    using_long_TOF = TOF > TOF_min;
    
    % implement fsolve
    f = @(a)Lamberts_root(a, mu, s, c, TOF, TOF_min, using_long_transfer_angle);
    delta_a = 10000;
    a0 = a_min + delta_a;
    % tolerance = 10^-12;
    tolerance = 10^-12;
%     options = optimoptions('fsolve','Display','iter', 'FunctionTolerance', tolerance, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
%     options = optimoptions('fsolve', 'FunctionTolerance', tolerance, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
    options = optimoptions('fsolve', 'Display','off', 'FunctionTolerance', tolerance, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
    [a,~,~,output] = fsolve(f, a0, options);
    iterations = output.iterations;

    [alpha, beta] = getAlphaBeta(c, s, a, using_long_TOF, using_long_transfer_angle);
    p = (4*a*(s-R1)*(s-R2)/(c^2))*((sin((alpha+beta)/2))^2);

    % added to output v1 and v2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % e = sqrt(1 - p/a); % e >= 0 -> no sign check needed
    % trueAnomaly1 = acos(1/e*(p/R1 - 1));
    % trueAnomaly2 = acos(1/e*(p/R2 - 1));
    transfer_angle = acos(dot(r1_N,r2_N)/R1/R2);
    if using_long_transfer_angle
        transfer_angle = 2*pi - transfer_angle;
    end
    % transfer_angle_deg = transfer_angle*180/pi;
    % tol = 10^-12;
    % if abs(trueAnomaly2 - trueAnomaly1 - transfer_angle) > tol
    %     if abs(trueAnomaly2 - (2*pi-trueAnomaly1) - transfer_angle) < tol
    %         trueAnomaly1 = 2*pi-trueAnomaly1;
    %     elseif abs(2*pi-trueAnomaly2 - trueAnomaly1 - transfer_angle) < tol
    %         trueAnomaly2 = 2*pi-trueAnomaly2;
    %     elseif abs(2*pi-trueAnomaly2 - (2*pi-trueAnomaly1) - transfer_angle) < tol
    %         trueAnomaly2 = -trueAnomaly2;
    %         trueAnomaly1 = -trueAnomaly1;
    %     else
    %         fprintf("True anomalies do not add to transfer angle")
    %     end
    % end
    % trueAnomaly1_deg = trueAnomaly1*180/pi;
    % trueAnomaly2_deg = trueAnomaly2*180/pi;
    
    f = 1 - R2/p*(1-cos(transfer_angle));
    g = R2*R1/sqrt(mu*p)*sin(transfer_angle);
    fdot = sqrt(mu/p)*tan(transfer_angle/2)*((1-cos(transfer_angle))/p - 1/R2 - 1/R1);
    gdot = 1-(R1/p)*(1 - cos(transfer_angle));
    v1_N = (r2_N - f*r1_N)/g;
    v2_N = fdot*r1_N + gdot*v1_N;
end

function f = Lamberts_root(a, mu, s, c, TOF, TOF_min, using_long_transfer_angle)
    n = sqrt(mu/a^3);
    using_long_TOF = TOF > TOF_min;
    [alpha, beta] = getAlphaBeta(c, s, a, using_long_TOF, using_long_transfer_angle);
    f = 1/n*(alpha-beta-(sin(alpha)-sin(beta))) - TOF;
end

function [c,s] = getChordAndSemiPerimeter(R1_SCI, R2_SCI, using_long_transfer_angle)
    R1 = norm(R1_SCI);
    R2 = norm(R2_SCI);
    delta_theta = acos(dot(R1_SCI,R2_SCI)/R1/R2); % want transfer angle less than 180
    if using_long_transfer_angle
        delta_theta = 2*pi - delta_theta;
    end
    c = sqrt(R1^2 + R2^2 - 2*R1*R2*cos(delta_theta));
    s = 1/2*(R1 + R2 + c);
end

function [alpha, beta] = getAlphaBeta(c, s, a, using_long_TOF, using_long_transfer_angle)
    alpha = 2*asin(sqrt(s/2/a));
    beta = 2*asin(sqrt((s-c)/2/a));
    alpha = real(alpha);
    beta = real(beta);
    if using_long_TOF
        alpha = 2*pi-alpha;
    end
    if using_long_transfer_angle
       beta = -beta; 
    end
end

