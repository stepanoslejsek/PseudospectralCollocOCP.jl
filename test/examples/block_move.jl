using PseudospectralCollocOCP

function dynamics(x, u, t)
    return [x[2], u[1]]
end

# Using tanh smoothing for abs value
α = 0.01
running_cost = (x, u, t) -> u[1] * x[2] * tanh(u[1] * x[2] / α)

x0 = [0.0, 0.0]
xf = [1.0, 0.0]

umax = 10.0
control_bounds = [(-umax, umax)]
tf = 1.0

t, x, u, obj = solve_ocp(N=150, tf=tf, x0=x0, xf=xf, dynamics=dynamics, running_cost=running_cost, nx=2, nu=1, control_bounds=control_bounds)

fig = Figure()
ax1 = Axis(fig[1, 1], ylabel="x [m]")
ax2 = Axis(fig[2, 1], ylabel="v [m/s]")
ax3 = Axis(fig[3, 1], xlabel="t [s]", ylabel="F [N]")
scatterlines!(ax1, t, x[1, :], color=:blue)
scatterlines!(ax2, t, x[2, :], color=:blue)
scatterlines!(ax3, t, u[1, :], color=:blue)
fig
