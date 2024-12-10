using OrdinaryDiffEq, LinearAlgebra, Plots, Statistics

function rv2orbitEls(posVel::Vector{Float64}, mu::Float64, using_trueAnomaly::Int = 0)
    # Canonical basis vectors
    X_N = [1.0, 0.0, 0.0]
    Y_N = [0.0, 1.0, 0.0]
    Z_N = [0.0, 0.0, 1.0]

    # Extract position and velocity
    r_N = posVel[1:3]
    v_N = posVel[4:6]
    r = norm(r_N)

    # Angular momentum
    h_N = cross(r_N, v_N)
    h = norm(h_N)
    p = h^2 / mu

    # Eccentricity vector and magnitude
    e_N = cross(v_N, h_N) / mu - r_N / r
    e = norm(e_N)
    a = p / (1 + e^2)

    # Line of nodes
    n_N = cross(Z_N, h_N)
    n = norm(n_N)
    OMEGA = acos(dot(n_N, X_N) / n)
    if dot(n_N, Y_N) < 0
        OMEGA = -OMEGA
    end

    omega = acos(dot(n_N, e_N) / (n * e))
    if dot(e_N, Z_N) < 0
        omega = -omega
    end

    i = acos(dot(h_N, Z_N) / h)

    nu = acos(dot(e_N, r_N) / (e * r))
    if dot(r_N, v_N) < 0
        nu = -nu
    end
    nu = real(nu)  # Ensure it's real

    if using_trueAnomaly == 1
        orbitEls = vcat(a, e, i, OMEGA, omega, nu)
    else
        M = atan2(-sqrt(1 - e^2) * sin(nu), -e - cos(nu)) +
            pi - e * sqrt(1 - e^2) * sin(nu) / (1 + e * cos(nu))
        orbitEls = vcat(a, e, i, OMEGA, omega, M)
    end

    theta = nu + omega
    ON = R3(theta) * R1(i) * R3(OMEGA)
    NO = transpose(ON)

    return orbitEls, NO
end

# Helper functions for rotation matrices
function R1(angle::Float64)
    return [
        1.0 0.0 0.0;
        0.0 cos(angle) -sin(angle);
        0.0 sin(angle) cos(angle)
    ]
end

function R3(angle::Float64)
    return [
        cos(angle) -sin(angle) 0.0;
        sin(angle) cos(angle) 0.0;
        0.0 0.0 1.0
    ]
end

# Canonical units from Earth's surface
TU = 806.812  # s
DU = 6378.140  # km
mu_canonical = 1.0  # Canonical units
mu = 3.986004415e5  # km^3/s^2

Fmax = 0.01  # Canonical units (9.8e-5 km/s^2)

# Initial conditions
r0 = [-0.70545852988580, -0.73885031681775, -0.40116299069586]
v0 = [0.73122658145185, -0.53921753373056, -0.29277123328399]
x0 = vcat(r0, v0)
orbitEls_t0,_ = rv2orbitEls(x0, mu_canonical, 1);

# Target angular momentum and Laplace vector
L_targ = [0.0, 0.0, 2.56612389857378]
A_targ = [0.0, 0.0, 0.0]
orbitEls_targ = [42000/DU; 0; 0; 0; 0; 0];

# Tuning parameters
epsilon = 0.00001
k = 2.0

# Julia implementation of control_F function
function control_F(r::Array{Float64, 2}, v::Array{Float64, 2}, L_targ::Array{Float64, 1}, A_targ::Array{Float64, 1}, Fmax::Float64, k::Float64, epsilon::Float64, mu::Float64)
    m, n = size(r)
    F = Array{Float64, 2}(undef, m, n)
    for i in 1:n
        L = cross(r[:, i], v[:, i])
        A = cross(v[:, i], L) - mu * r[:, i] / norm(r[:, i])
        delta_L = L - L_targ
        delta_A = A - A_targ
        G = -(cross(k * delta_L, r[:, i]) + cross(L, delta_A) + cross(cross(delta_A, v[:, i]), r[:, i]))
        if norm(G) < epsilon * Fmax
            F[:, i] = G / epsilon
        else
            F[:, i] = Fmax * G / norm(G)
        end
    end
    return F
end

# Julia implementation of lowThrustDynamics function
function lowThrustDynamics(X::Array{Float64, 1}, L_targ::Array{Float64, 1}, A_targ::Array{Float64, 1}, Fmax::Float64, k::Float64, epsilon::Float64, mu::Float64)
    r = X[1:3]  # Extract position from state vector
    v = X[4:6]  # Extract velocity from state vector
    dr = v      # Time derivative of position is velocity

    # Compute control force using control_F
    F = control_F(reshape(r, 3, 1), reshape(v, 3, 1), L_targ, A_targ, Fmax, k, epsilon, mu)[:, 1]
    dv = -mu / (norm(r)^3) * r + F  # Inertial acceleration

    return vcat(dr, dv)  # Concatenate derivatives for state vector
end

# Define dynamics
function dynamics!(du, u, p, t)
    L_targ, A_targ, Fmax, k, epsilon, mu = p
    du[:] = lowThrustDynamics(u, L_targ, A_targ, Fmax, k, epsilon, mu)
end

# Time span and solve
tspan = (0.0, 20 / TU * 3600)  # Canonical time
params = (L_targ, A_targ, Fmax, k, epsilon, mu_canonical)
prob = ODEProblem(dynamics!, x0, tspan, params)
sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12) #, maxiters=1e7)
# sol = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12)
# sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
# sol = solve(prob, DP5(), reltol=1e-12, abstol=1e-12)
# sol = solve(prob, AutoVern7(Rodas5()), reltol=1e-12, abstol=1e-12)


# Compute control forces
r_sol = hcat([sol.u[i][1:3] for i in 1:length(sol)]...)
v_sol = hcat([sol.u[i][4:6] for i in 1:length(sol)]...)
F_sol = [control_F(r_sol[:, i:i], v_sol[:, i:i], L_targ, A_targ, Fmax, k, epsilon, mu_canonical)[:, 1] for i in 1:size(r_sol, 2)]
F_sol = hcat(F_sol...)

# Convert to physical units
t_physical = sol.t .* TU
F_physical = F_sol .* DU / (TU^2)

# gr(size=(600,600))

function plot_iso3d(xs, ys, zs; lw=3, ls=:solid, label=false)
    # condition data for nearly isometric 3D plot 
    x12, y12, z12 = extrema(xs), extrema(ys), extrema(zs)
    d = maximum([diff([x12...]),diff([y12...]),diff([z12...])])[1] / 2
    xm, ym, zm = mean(x12),  mean(y12),  mean(z12) 

    # plot data
    Plots.plot!(; xlabel="x",ylabel="y",zlabel="z", aspect_ratio=:equal, grid=:true, size=(600,600))
    Plots.plot!(xlims=(xm-d,xm+d), ylims=(ym-d,ym+d), zlims=(zm-d,zm+d))
    Plots.plot!(xs, ys, zs,lw=lw,ls=ls,label=label)
end

function plotOrbit(num_points::Int, orbitEls::Vector{Float64}, lineStyle::Symbol, lineWidth::Float64, labelStr::String)
    # Unpack orbital elements
    a, e, i, OMEGA, omega, _ = orbitEls  # Assuming orbitEls is in the correct order

    # True anomaly range
    nu_range = range(0, 2π, length=num_points)

    # Preallocate arrays for Cartesian coordinates
    x = zeros(Float64, num_points)
    y = zeros(Float64, num_points)
    z = zeros(Float64, num_points)

    for j in 1:num_points
        # Current true anomaly
        nu_current = nu_range[j]

        # Radius at current true anomaly
        r = (a * (1 - e^2)) / (1 + e * cos(nu_current))

        # Position in the orbital plane
        x_orbital = r * cos(nu_current)
        y_orbital = r * sin(nu_current)

        # Convert to 3D Cartesian coordinates
        x[j] = (cos(omega) * cos(OMEGA) - sin(omega) * sin(OMEGA) * cos(i)) * x_orbital +
               (-sin(omega) * cos(OMEGA) - cos(omega) * sin(OMEGA) * cos(i)) * y_orbital
        y[j] = (cos(omega) * sin(OMEGA) + sin(omega) * cos(OMEGA) * cos(i)) * x_orbital +
               (-sin(omega) * sin(OMEGA) + cos(omega) * cos(OMEGA) * cos(i)) * y_orbital
        z[j] = (sin(omega) * sin(i)) * x_orbital +
               (cos(omega) * sin(i)) * y_orbital
    end

    # Plot the orbit in 3D
    # plot3d!(x, y, z, label=labelStr, linestyle=lineStyle, linewidth=lineWidth, aspect_ratio=:equal)
    plot_iso3d(x, y, z; lw=lineWidth, ls=lineStyle, label=labelStr)
end

# gr(size=(600,600))
# Plot trajectory
# plot3d(r_sol[1, :], r_sol[2, :], r_sol[3, :], label="Trajectory", linewidth=1.5, aspect_ratio=:equal)
p = plot()
plot_iso3d(r_sol[1, :], r_sol[2, :], r_sol[3, :]; lw=1.5, label="Trajectory")
# plot3d(r_sol[:,1], r_sol[:,2], r_sol[:,3], label="Trajectory", linewidth=1.5, aspect_ratio=1.0)
plotOrbit(100, orbitEls_t0, :dash, 1.5, "Initial Orbit")
plotOrbit(100, orbitEls_targ, :dash, 1.5, "Target Orbit")

# Convert time to hours and force to m/s^2
time_hours = t_physical / 3600  # Convert seconds to hours
F_mps2 = F_physical * 1e3       # Convert from km/s^2 to m/s^2

# # Plot each component of the control force
# plot(time_hours, F_mps2[1, :], label="Fx", linewidth=2)
# plot!(time_hours, F_mps2[2, :], label="Fy", linewidth=2)
# plot!(time_hours, F_mps2[3, :], label="Fz", linewidth=2)

# Add labels and legend
# xlabel!("time (hours)")
# ylabel!("Control Force (m/s^2)")