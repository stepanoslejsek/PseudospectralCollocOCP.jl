using PseudospectralCollocOCP

m1 = 1.0
m2 = 0.3
l = 0.5
d = 1.0
g = 9.81
tf = 2.0
umax = 20.0
dmax = 2.0

function dynamics(x, u, t)
    # [x, theta, dot_x, dot_theta]
    ẋ = x[3]
    θ̇ = x[4]
    ẍ = (l * m2 * sin(x[2]) * x[4]^2 + u[1] + m2 * g * cos(x[2]) * sin(x[2])) / (m1 + m2 * (1 - cos(x[2])^2))
    θ̈ = -(l * m2 * cos(x[2]) * sin(x[2]) * x[4]^2 + u[1] * cos(x[2]) + (m1 + m2) * g * sin(x[2])) / (l * m1 + l * m2 * (1 - cos(x[2])^2))
    return [ẋ, θ̇, ẍ, θ̈]
end

running_cost = (x, u, t) -> u[1]^2

nx = 4
nu = 1
x0 = zeros(nx)
xf = [d, pi, 0, 0]

control_bounds = [(-umax, umax)]
state_bounds = [(-dmax, dmax), nothing, nothing, nothing]

t, x, u, obj = solve_ocp(N=50, tf=tf, x0=x0, xf=xf, dynamics=dynamics, running_cost=running_cost, nx=nx, nu=nu)

fig = Figure()
ax1 = Axis(fig[1, 1], ylabel="x [m]")
ax2 = Axis(fig[2, 1], ylabel="theta [rad]")
ax3 = Axis(fig[3, 1], xlabel="t [s]", ylabel="F [N]")
scatterlines!(ax1, t, x[1, :], color=:blue)
scatterlines!(ax2, t, x[2, :], color=:blue)
scatterlines!(ax3, t, u[1, :], color=:blue)
fig
