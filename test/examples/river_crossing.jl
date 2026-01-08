using PseudospectralCollocOCP
using GLMakie

Rmax = 3.0
wmax = 1.0
R(w) = Rmax * exp(-(w / wmax)^2)
U(x, y) = R(y - sin(x)) / sqrt(1 + cos(x)^2)
V(x, y) = R(y - sin(x)) * cos(x) / sqrt(1 + cos(x)^2)

function dynamics(x, u, t)
    dot_x = x[3] * cos(x[4]) + U(x[1], x[2])
    dot_y = x[3] * sin(x[4]) + V(x[1], x[2])
    dot_v = u[1]
    dot_theta = u[2]
    return [dot_x, dot_y, dot_v, dot_theta]
end

# Time-optimal control
running_cost = (x, u, t) -> 1

x0 = [0.0, -1.0, nothing, nothing]
xf = [2 * pi, 1.0, nothing, nothing]

southern_bank = (x, u, t) -> sin.(x[1, :]) .- 1.0
northern_bank = (x, u, t) -> sin.(x[1, :]) .+ 1.0

state_bounds = [nothing, (southern_bank, northern_bank), (0.0, 4.0), nothing]
control_bounds = [(-100.0, 100.0), (-5.0, 5.0)]

t, x, u, obj, x_interp, u_interp = solve_ocp(N=50, tf=nothing, x0=x0, xf=xf, dynamics=dynamics, running_cost=running_cost, nx=4, nu=2, state_bounds=state_bounds, control_bounds=control_bounds)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="x [m]", ylabel="y [m]")
scatterlines!(ax, x[1, :], x[2, :], color=:blue)
lines!(ax, x[1, :], southern_bank(x, 0, 0), color=:black, linewidth=3)
lines!(ax, x[1, :], northern_bank(x, 0, 0), color=:black, linewidth=3)

N = 500
t_dense = collect(range(t[1], t[end], N))
x_dense = hcat(x_interp.(t_dense)...)
scatterlines!(ax, x_dense[1, :], x_dense[2, :], color=:red)
fig
