using FastGaussQuadrature

function legendre_polynomial(n::Int, x::Float64)
    # Bonnet's recursion formula
    if n == 0
        return 1.0
    elseif n == 1
        return x
    else
        P0 = 1.0
        P1 = x
        for k in 2:n
            P2 = ((2k - 1) * x * P1 - (k - 1) * P0) / k
            P0 = P1
            P1 = P2
        end
        return P1
    end
end

function LGL_collocation(N::Int)
    nodes, weights = gausslobatto(N)

    D = zeros(N, N)

    # Differentiation matrix
    for i in 1:N
        for j in 1:N
            if i != j
                Pn_xi = legendre_polynomial(N - 1, nodes[i])
                Pn_xj = legendre_polynomial(N - 1, nodes[j])
                D[i, j] = Pn_xi / (Pn_xj * (nodes[i] - nodes[j]))
            end
        end
    end
    for i in 1:N
        D[i, i] = -sum(D[i, setdiff(1:N, i)])
    end
    return nodes, weights, D
end
