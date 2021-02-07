include("Graph.jl")

using LinearAlgebra
using SparseArrays

function Degree(G)
    C = zeros(G.n)
    for (ID, u, v, w) in G.E
        C[u] += w
        C[v] += w
    end
    return C
end

function spA(G)
    Is = zeros(Int, G.m*2)
    Js = zeros(Int, G.m*2)
    Vs = zeros(G.m*2)
    for (ID, u, v, w) in G.E
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = w
        Is[G.m + ID] = v
        Js[G.m + ID] = u
        Vs[G.m + ID] = w
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end

function alpP(G, alp)
    d = Degree(G)
    Is = zeros(Int, G.m*2)
    Js = zeros(Int, G.m*2)
    Vs = zeros(G.m*2)
    for (ID, u, v, w) in G.E
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = alp * w / d[v]
        Is[G.m + ID] = v
        Js[G.m + ID] = u
        Vs[G.m + ID] = alp * w / d[u]
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end

function BetweennessCentrality(G)
    g = getAdjacentList(G)
    C = zeros(G.n)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    Q = zeros(Int32, G.n+10)
    delta = zeros(G.n)
    for s = 1 : G.n
        foreach(i -> p[i] = [], 1 : G.n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                if w != s
                    C[w] += delta[w]
                end
            end
        end

    end

    return C
end

function selectPivots(G, eps)
    deg = Degree(G)
    nl = zeros(Int, G.n)
    foreach(i -> nl[i] = i, 1 : G.n)
    sort!(nl, by = x->-deg[x]);
    pivots = []
    kk = round(Int, log2(G.n)/eps^2)
    foreach(i -> push!(pivots, nl[i]), 1 : kk)
    return pivots
end

function ApproxBetweennessCentrality(G; eps = 0.5)
    g = getAdjacentList(G)
    pivots = selectPivots(G, eps)
    C = zeros(G.n)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    Q = zeros(Int32, G.n+10)
    delta = zeros(G.n)
    for s in pivots
        foreach(i -> p[i] = [], 1 : G.n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                if w != s
                    C[w] += delta[w]
                end
            end
        end

    end

    return C
end

function ClosenessCentrality(G)
    g = getAdjacentList(G)
    C = zeros(G.n)
    d = zeros(Int32, G.n)
    Q = zeros(Int32, G.n+10)
    for s = 1 : G.n
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
            end
        end

        C[s] = sum(d)
    end

    foreach(i -> C[i] = 1.0 / C[i], 1 : G.n)

    return C
end

function ForestDistanceClosenessCentrality(G)
    # get Forest Matrix
    IpL = zeros(G.n, G.n)
    for (ID, u, v, w) in G.E
        IpL[u, u] += w
        IpL[v, v] += w
        IpL[u, v] -= w
        IpL[v, u] -= w
    end
    foreach(i -> IpL[i, i] += 1, 1 : G.n)
    W = inv(IpL)

    # calculate FDC
    trace = tr(W)
    C = zeros(G.n)
    foreach(i -> C[i] = G.n / (G.n * W[i, i] + trace - 2), 1 : G.n)
    #foreach(i -> C[i] = W[i, i], 1 : G.n)

    return C
end

function InformationCentrality(G)
    # get pseudo inverse of L
    L = zeros(G.n, G.n)
    for (ID, u, v, w) in G.E
        L[u, u] += w
        L[v, v] += w
        L[u, v] -= w
        L[v, u] -= w
    end
    L .+= (1 / G.n)
    Lp = inv(L)
    Lp .-= (1 / G.n)

    # calculate Information Centrality
    trace = tr(Lp)
    C = zeros(G.n)
    foreach(i -> C[i] = G.n / (G.n * Lp[i, i] + trace), 1 : G.n)
    #foreach(i -> C[i] = Lp[i, i], 1 : G.n)

    return C
end

function PageRank(G; alpha = 0.85)
    aP = alpP(G, alpha)
    C = zeros(G.n)
    C[1] = 1.0
    adC = ((1.0 - alpha) / G.n) * ones(G.n)

    while true
        pC = copy(C)
        C = aP * C + adC
        if norm(C - pC) < 1e-12
            break
        end
    end

    return C
end

function LeaderRank(G)
    n = G.n + 1
    g = Array{Array{Int32, 1}, 1}(undef, n)
    foreach(i -> g[i] = [], 1 : n)
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    for i = 1 : n-1
        push!(g[n], i)
        push!(g[i], n)
    end
    d = zeros(Int32, n)
    foreach(i -> d[i] = size(g[i], 1), 1 : n)

    C = zeros(n)
    err = zeros(n)
    pC = zeros(n)
    C[1] = 1.0

    while true
        foreach(i -> pC[i] = C[i], 1 : n)
        for u = 1 : n
            C[u] = 0
            for v in g[u]
                C[u] += (pC[v] / d[v])
            end
        end
        foreach(i -> err[i] = abs(C[i] - pC[i]), 1 : n)
        maxErr = err[argmax(err)]
        if maxErr < 1e-9
            break
        end
    end

    ansC = zeros(G.n)
    foreach(i -> ansC[i] = C[i], 1 : G.n)

    return ansC
end

function LeverageCentrality(G)
    d = Degree(G)
    C = zeros(G.n)
    for (ID, u, v, w) in G.E
        C[u] += ((d[u] - d[v]) / (d[u] + d[v]))
        C[v] += ((d[v] - d[u]) / (d[u] + d[v]))
    end
    foreach(i -> C[i] /= d[i], 1 : G.n)
    return C
end

function EigenvectorCentrality(G)
    C = zeros(G.n)
    C[1] = 1.0
    A = spA(G)

    while true
        pC = copy(C)
        C = A * C
        C ./= C[argmax(C)]
        if norm(C - pC) < 1e-9
            break
        end
    end

    return C
end

function KatzCentrality(G; fac = 0.85)
    A = spA(G)
    getmax(X) = X[argmax(X)]
    x = ones(G.n)
    xc = ones(G.n)
    tt = 0.0
    lmd = 0.0
    while true
        xc = copy(x)
        x = A * x
        tmp = getmax(x)
        x ./= tmp
        if abs(tmp-tt)<1e-6
            lmd = tmp
            break
        end
        tt = tmp
    end
    alp = fac / lmd
    C = zeros(G.n)
    addC = ones(G.n)
    while getmax(addC) > 1e-9
        addC = A * addC
        addC .*= alp
        C = C + addC
    end

    return C
end

function SpanningTreeCentrality(G)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    d = zeros(Int32, G.n)
    foreach(i -> d[i] = size(g[i], 1), 1 : G.n)

    L = zeros(G.n-1, G.n-1)
    for (ID, u, v, w) in G.E
        if u < G.n
            L[u, u] += w
        end
        if v < G.n
            L[v, v] += w
        end
        if (u < G.n) && (v < G.n)
            L[u, v] -= w
            L[v, u] -= w
        end
    end
    nL = inv(L)

    I(N) = Diagonal(ones(N))

    C = zeros(G.n)

    for u = 1 : G.n
        if d[u] == 1
            C[u] = 0.0
        else
            U = spzeros(d[u]-1, G.n-1)
            for i = 1 : d[u]-1
                v = g[u][i]
                if u == G.n
                    U[i, v] = 1
                elseif v == G.n
                    U[i, u] = 1
                else
                    U[i, u] = 1
                    U[i, v] = -1
                end
            end
            C[u] = 1.0 - d[u] * det(I(d[u]-1) - U * nL * U')
        end
    end

    return C
end

function STnumber(G)
    L = zeros(G.n-1, G.n-1)
    for (ID, u, v, w) in G.E
        L[u, u] += w
        if v < G.n
            L[v, v] += w
            L[u, v] -= w
            L[v, u] -= w
        end
    end
    return det(L)
end

function reduced(G, sv)
    V = zeros(Int32, G.n-1)
    ori = zeros(Int32, G.n)
    for i = 1 : sv-1
        V[i] = i
        ori[i] = i
    end
    for i = sv+1 : G.n
        V[i-1] = i
        ori[i] = i-1
    end
    edge = Set{Tuple{Int32, Int32}}()
    for (ID, u, v, w) in G.E
        if (u != sv) && (v != sv)
            push!(edge, (ori[u], ori[v]))
        end
    end
    m = length(edge)
    E = Array{Tuple{Int32, Int32, Int32, Float64}, 1}(undef, m)

    ID = 0
    n = G.n-1
    for (u, v) in edge
        ID = ID + 1
        E[ID] = (ID, u, v, 1.0)
    end

    return Graph(n, m, V, E)
end

function STCdp(G)
    C = zeros(G.n)
    all = STnumber(G)
    d = Degree(G)
    for i = 1 : G.n
        C[i] = (all - d[i] * STnumber(reduced(G, i))) / all
    end
    return C
end
