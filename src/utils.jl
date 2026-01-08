using FastGaussQuadrature

t2τ(t, t0, tf) = 2 ./ (tf .- t0) .* t .- (tf .+ t0) ./ (tf .- t0)
τ2t(τ, t0, tf) = (tf .+ t0) ./ 2 .+ (tf .- t0) ./ 2 .* τ

function legendre_polynomial(n::Int, x::Float64)
    # Returns P_n(x) and P_n'(x)
    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return x, 1.0
    else
        P0 = 1.0   # P₀
        P1 = x     # P₁

        for k in 2:n
            P2 = ((2k - 1) * x * P1 - (k - 1) * P0) / k
            P0 = P1
            P1 = P2
        end

        Pn = P1
        Pn_minus1 = P0

        # Derivative using exact identity
        if abs(abs(x) - 1.0) < 1e-14
            dPn = 0.5 * n * (n + 1) * (x == 1.0 ? 1.0 : (-1.0)^(n + 1))
        else
            dPn = n * (x * Pn - Pn_minus1) / (x^2 - 1.0)
        end

        return Pn, dPn
    end
end

function LGL_collocation(N::Int)
    nodes, weights = gausslobatto(N)

    D = zeros(N, N)

    # Differentiation matrix
    for i in 1:N
        for j in 1:N
            if i != j
                Pn_xi = legendre_polynomial(N - 1, nodes[i])[1]
                Pn_xj = legendre_polynomial(N - 1, nodes[j])[1]
                D[i, j] = Pn_xi / (Pn_xj * (nodes[i] - nodes[j]))
            end
        end
    end
    for i in 1:N
        D[i, i] = -sum(D[i, setdiff(1:N, i)])
    end
    return nodes, weights, D
end

function compute_barycentric_weights(nodes::Vector{Float64}, N::Int)
    # For LGL points, barycentric weights can be computed from Legendre polynomials
    w = zeros(length(nodes))
    for k in 1:length(nodes)
        Pn = legendre_polynomial(N - 1, nodes[k])[1]
        w[k] = 1.0 / Pn
    end
    return w
end

function lagrange_polynomial(j::Int, t::Float64, nodes::Vector{Float64}, weights_bary::Vector{Float64})
    N = length(nodes)

    # Handle the case where t coincides with a node
    for i in 1:N
        if abs(t - nodes[i]) < 1e-14
            return i == j ? 1.0 : 0.0
        end
    end

    # Barycentric formula
    numerator = weights_bary[j] / (t - nodes[j])
    denominator = sum(weights_bary[k] / (t - nodes[k]) for k in 1:N)

    return numerator / denominator
end
