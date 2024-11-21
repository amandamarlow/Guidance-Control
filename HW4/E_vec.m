function [E] = E_vec(x_tf_N, rp_N, vp_N, mu)

    r_tf_N = x_tf_N(1:3);
    r_tf = norm(r_tf_N);
    v_tf_N = x_tf_N(4:6);

    rp = norm(rp_N);
    vp = norm(vp_N);
    rp_hat = rp_N/rp;
    vp_hat = vp_N/vp;
    h_hat = cross(rp_hat, vp_hat);
    p = vp^2*rp^2/mu;
    e_N = cross(v_N,h_N)/mu - r_N/r;

    nu = atan2(dot(r_tf_N,vp_hat), dot(r_tf_N,rp_hat));
    rR = p/(1+e*cos(nu));
    vR_N = -1/rp*sqrt(mu/p)*sin(nu)*rp_N + (1-rp/p*(1-cos(nu)))*vp_N;
    E(1:3) = vR_N - v_tf_N;
    E(4) = rR - r_tf;
    E(5) = dot(h_hat, r_tf_N);
    g = -mu/r_tf^3*r_tf_N;
    lambda_tf = A_N + B_N*tf;
    E(6) = dot(lambda_tf, g) - dot(B_N, v_tf_N);
    E(7) = sqrt(lambda_tf'*lambda_tf)-1;

end

