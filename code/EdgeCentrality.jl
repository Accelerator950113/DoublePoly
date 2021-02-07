include("Graph.jl")
include("NodeCentrality.jl")

using SparseArrays
using LinearAlgebra

function spA(G)
    Is = zeros(Int, 2*G.m)
    Js = zeros(Int, 2*G.m)
    Vs = zeros(2*G.m)
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

function EdgeKatz(G; fac1 = 0.85)
    KC = KatzCentrality(G, fac=fac1)
    EK = zeros(G.m)
    for (ID, u, v, w) in G.E
        EK[ID] = (1.0 + KC[u]) * (1.0 + KC[v])
    end
    return EK
end

function getLp(G)
    L = zeros(G.n, G.n)
    for (ID, u, v, w) in G.E
        L[u, v] -= w
        L[v, u] -= w
        L[u, u] += w
        L[v, v] += w
    end
    L .+= (1.0 / G.n)
    Lp = inv(L)
    Lp .-= (1.0 / G.n)
    return Lp
end

function BiharmonicDistanceRelatedCentrality(G)
    Lp = getLp(G)
    C = zeros(G.m)
    for (ID, u, v, w) in G.E
        C[ID] = norm(Lp[:, u] - Lp[:, v])^2
    end
    return C
end

function SpanningEdgeCentrality(G)
    Lp = getLp(G)
    C = zeros(G.m)
    for (ID, u, v, w) in G.E
        C[ID] = Lp[u, u] + Lp[v, v] - 2*Lp[u, v]
    end
    return C
end

function EdgeBetweenness(G)
    g = getAdjacentList(G)
    idx = Dict{Tuple{Int32, Int32}, Int32}()
    for (ID, u, v, w) in G.E
        idx[(u, v)] = ID
        idx[(v, u)] = ID
    end
    Cb = zeros(G.m)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    d = zeros(Int32, G.n)
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
                Cb[idx[(v, w)]] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
            end
        end

    end

    return Cb
end

function ApproxEdgeBetweenness(G; eps = 0.5)
    g = getAdjacentList(G)
    idx = Dict{Tuple{Int32, Int32}, Int32}()
    for (ID, u, v, w) in G.E
        idx[(u, v)] = ID
        idx[(v, u)] = ID
    end
    pivots = selectPivots(G, eps)
    Cb = zeros(G.m)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    d = zeros(Int32, G.n)
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
                Cb[idx[(v, w)]] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
            end
        end

    end

    return Cb
end

function getW(G)
    IpL = zeros(G.n, G.n)
    for (ID, u, v, w) in G.E
        IpL[u, v] -= w
        IpL[v, u] -= w
        IpL[u, u] += w
        IpL[v, v] += w
    end
    foreach(i -> IpL[i, i] += 1, 1 : G.n)
    return inv(IpL)
end

function ForestCentrality(G)
    W = getW(G)
    C = zeros(G.m)
    for (ID, u, v, w) in G.E
        C[ID] = (W[u, u] + W[v, v] - 2*W[u, v]) / W[u, v]
    end
    return C
end
