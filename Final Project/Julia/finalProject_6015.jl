using DifferentialEquations, LinearAlgebra, Plots

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

# Target angular momentum and Laplace vector
L_targ = [0.0, 0.0, 2.56612389857378]
A_targ = [0.0, 0.0, 0.0]

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
sol = solve(prob, reltol=1e-12, abstol=1e-12)

# Compute control forces
r_sol = hcat([sol.u[i][1:3] for i in 1:length(sol)]...)
v_sol = hcat([sol.u[i][4:6] for i in 1:length(sol)]...)
F_sol = [control_F(r_sol[:, i:i], v_sol[:, i:i], L_targ, A_targ, Fmax, k, epsilon, mu_canonical)[:, 1] for i in 1:size(r_sol, 2)]
F_sol = hcat(F_sol...)

# Convert to physical units
t_physical = sol.t .* TU
F_physical = F_sol .* DU / (TU^2)

# Plot trajectory
plot3d(r_sol[1, :], r_sol[2, :], r_sol[3, :], label="Trajectory", linewidth=1.5)
