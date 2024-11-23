function [E] = ode45wrapper(t, x0, lambdaAug, maxThrust, targetOrbitEls, mu)

    % A_N = lambdaAug(4:6);
    % B_N = lambdaAug(1:3);
    A = lambdaAug(4:6);
    B = lambdaAug(1:3);
    t_go = lambdaAug(end);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,x] = ode45(@(t,x) bilinearTangent_dynamics(t, x, mu, A, B, maxThrust), [t t+t_go], x0, options);

    x_tf_N = x(end,:)';
    r_tf_N = x_tf_N(1:3);
    r_tf = norm(r_tf_N);
    v_tf_N = x_tf_N(4:6);
    [~, NO_tf] = rv2orbitEls(x_tf_N, mu, 1);
    r_tf_Otf = NO_tf'*r_tf_N;
    v_tf_Otf = NO_tf'*v_tf_N;

    % rp = norm(rp_N);
    % vp = norm(vp_N);
    % rp_hat = rp_N/rp;
    % vp_hat = vp_N/vp;
    % h_hat = cross(rp_hat, vp_hat);
    % p = vp^2*rp^2/mu;
    % e_N = cross(v_N,h_N)/mu - r_N/r;
    
    a = targetOrbitEls(1);
    e = targetOrbitEls(2);
    rp = a*(1-e);
    vp = sqrt(2*mu/rp - mu/a);
    rp_O = [rp; 0; 0];
    vp_O = [0; vp; 0];
    [~, NO_periapsis] = orbitEls2xv(targetOrbitEls, mu);
    rp_N = NO_periapsis*rp_O;
    vp_N = NO_periapsis*vp_O;
    rp_hat = rp_N/rp;
    vp_hat = vp_N/vp;
    % h = sqrt(mu*a*(1-e^2));
    h_N = cross(rp_N, vp_N);
    h_hat = h_N/norm(h_N);
    p = vp^2*rp^2/mu;

    % OO_periapsis2tf = NO_tf'*NO_periapsis;
    % rp_Otf = OO_periapsis2tf*rp_O;
    % vp_Otf = OO_periapsis2tf*vp_O;
    % rp_hat = rp_Otf/rp;
    % vp_hat = vp_Otf/vp;
    % h = sqrt(mu*a*(1-e^2));
    % h_Otf = cross(rp_Otf, vp_Otf);
    % h_hat = h_Otf/h;
    
    E = NaN(7,1);
    nu = atan2(dot(r_tf_N,vp_hat), dot(r_tf_N,rp_hat));
    rR = p/(1+e*cos(nu));
    vR_N = -1/rp*sqrt(mu/p)*sin(nu)*rp_N + (1-rp/p*(1-cos(nu)))*vp_N;
    E(1:3) = vR_N - v_tf_N;
    E(4) = rR - r_tf;
    E(5) = dot(h_hat, r_tf_N);
    g = -mu/r_tf^3*r_tf_N;
    % lambda_tf = A_N + B_N*tf;
    % E(6) = dot(lambda_tf, g) - dot(B_N, v_tf_N);
    % E(7) = sqrt(lambda_tf'*lambda_tf)-1;
    lambda_tf = A + B*(t+t_go);
    % E(6) = dot(lambda_tf, g) - dot(NO_tf*B, v_tf_N);
    % E(6) = dot(lambda_tf, g) - dot(NO_tf*R3(-pi/2)*R1(-pi/2)*B, v_tf_N);
    E(6) = dot(lambda_tf, g) - dot(B, v_tf_N);
    E(7) = sqrt(lambda_tf'*lambda_tf)-1;

    % E = NaN(7,1);
    % nu = atan2(dot(r_tf_Otf,vp_hat), dot(r_tf_Otf,rp_hat));
    % rR = p/(1+e*cos(nu));
    % vR_Otf = -1/rp*sqrt(mu/p)*sin(nu)*rp_Otf + (1-rp/p*(1-cos(nu)))*vp_Otf;
    % E(1:3) = vR_Otf - v_tf_Otf;
    % E(4) = rR - r_tf;
    % E(5) = dot(h_hat, r_tf_Otf);
    % g = -mu/r_tf^3*r_tf_Otf;
    % % g = [-mu/r_tf^2, 0, 0];
    % lambda_tf = A + B*tf;
    % E(6) = dot(lambda_tf, g) - dot(B, v_tf_Otf);
    % E(7) = sqrt(dot(lambda_tf,lambda_tf)) - 1;

end

