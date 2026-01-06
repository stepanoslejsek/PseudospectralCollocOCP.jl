using JuMP
using Ipopt

function solve_ocp(;
    N::Int=50,
    t0::Float64=0.0,
    tf::Union{Nothing,Float64}=nothing,
    x0=nothing,
    xf=nothing,
    dynamics,
    running_cost=(x, u, t) -> 0.0,
    terminal_cost=(xf, tf) -> 0.0,
    state_bounds=nothing,
    control_bounds=nothing,
    nx::Int,
    nu::Int
)
    tau, w, D = LGL_collocation(N)

    model = Model(Ipopt.Optimizer)
    # Minimum time control
    if isnothing(tf)
        @variable(model, tf >= 0, start = 5.0)
    end

    t = (tf + t0) / 2 .+ (tf - t0) / 2 .* tau

    @variable(model, x[1:nx, 1:N])
    @variable(model, u[1:nu, 1:N])

    if !isnothing(state_bounds)
        for i in 1:nx
            if !isnothing(state_bounds[i])
                if typeof(state_bounds[i][1]) <: Number
                    # Static state bounds
                    set_lower_bound.(x[i, :], state_bounds[i][1])
                else
                    # Dynamic state bounds
                    @constraint(model, x[i, :] >= state_bounds[i][1](x, u, t))
                end
                if typeof(state_bounds[i][2]) <: Number
                    # Static state bounds
                    set_upper_bound.(x[i, :], state_bounds[i][2])
                else
                    # Dynamic state bounds
                    @constraint(model, x[i, :] <= state_bounds[i][2](x, u, t))
                end
            end
        end
    end

    if !isnothing(control_bounds)
        for i in 1:nu
            if !isnothing(control_bounds[i])
                if typeof(control_bounds[i][1]) <: Number
                    # Static control bounds
                    set_lower_bound.(u[i, :], control_bounds[i][1])
                else
                    # Dynamic control bounds
                    @constraint(model, u[i, :] >= control_bounds[i][1](x, u, t))
                end
                if typeof(control_bounds[i][2]) <: Number
                    # Static control bounds
                    set_upper_bound.(u[i, :], control_bounds[i][2])
                else
                    # Dynamic control bounds
                    @constraint(model, u[i, :] <= control_bounds[i][2](x, u, t))
                end
            end
        end
    end

    for i in 1:nx
        if !isnothing(x0[i])
            fix(x[i, 1], x0[i]; force=true)
        end
    end

    if !isnothing(xf)
        for i in 1:nx
            if !isnothing(xf[i])
                fix(x[i, N], xf[i]; force=true)
            end
        end
    end

    @objective(model, Min, (tf - t0) / 2 * sum(w[k] * running_cost(x[:, k], u[:, k], t[k]) for k in 1:N) + terminal_cost(x[:, N], tf))

    for k in 1:N
        for i in 1:nx
            @constraint(model, sum(D[k, j] * x[i, j] for j in 1:N) == (tf - t0) / 2 * dynamics(x[:, k], u[:, k], t[k])[i])
        end
    end

    optimize!(model)
    println("""
        termination_status = $(termination_status(model))
        primal_status      = $(primal_status(model))
        objective_value    = $(objective_value(model))
    """)
    return value.(t), value.(x), value.(u), objective_value(model)
end
