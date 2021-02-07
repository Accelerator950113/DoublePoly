include("Graph.jl")
include("NodeCentrality.jl")
include("EdgeCentrality.jl")

function NormalizedBetweennessCentrality(G)
    NC = ((G.n > 10000) || (G.m > 10000)) ? ApproxBetweennessCentrality(G) : BetweennessCentrality(G)
    EC = ((G.n > 10000) || (G.m > 10000)) ? ApproxEdgeBetweenness(G) : EdgeBetweenness(G)
    NC ./= sum(NC)
    EC ./= sum(EC)
    return NC, EC
end

function selectNodes(G, selectedNodes, selectedEdges, selectMode)
    NC, EC = NormalizedBetweennessCentrality(G)
    nl = zeros(Int, G.n)
    el = zeros(Int, G.m)
    foreach(i -> nl[i] = i, 1 : G.n)
    foreach(i -> el[i] = i, 1 : G.m)
    sort!(nl, by = x->-NC[x])
    sort!(el, by = x->-EC[x])

    Selected = zeros(Bool, G.n)
    if ('n' in selectMode)
        foreach(i -> Selected[nl[i]] = true, 1 : min(selectedNodes, G.n))
    end
    if ('e' in selectMode)
        for i = 1 : min(selectedEdges, G.m)
            Selected[G.E[el[i]][2]] = true
            Selected[G.E[el[i]][3]] = true
        end
    end

    return NC, EC, Selected
end

function SearchWay(g, status)
    n = size(g, 1)
    h = copy(status)
    Q = zeros(Int, n)
    front = 1
    rear = 0

    pushQ(x) = begin
        rear += 1
        Q[rear] = x
    end
    popQ() = begin
        tmp = Q[front]
        front += 1
        return tmp
    end

    for i = 1 : n
        (h[i] == 1) ? pushQ(i) : nothing
    end

    prior = zeros(Int, n)

    while front <= rear
        u = popQ()
        for v in g[u]
            if h[v] == 0
                h[v] = 1
                prior[v] = u
                pushQ(v)
            elseif h[v] == 2
                nodeList = Array{Int, 1}()
                push!(nodeList, v)
                push!(nodeList, u)
                p = prior[u]
                while p != 0
                    push!(nodeList, p)
                    p = prior[p]
                end
                return nodeList
            end
        end
    end

    return nothing
end

function KeepOne(G, Selected, NC, EC)
    g = getAdjacentList(G)
    F = zeros(Int, G.n)
    foreach(i -> F[i] = i, 1 : G.n)

    find(x) = begin
        if F[x] != x
            F[x] = find(F[x])
        end
        return F[x]
    end
    SetUnion(x, y) = F[find(x)] = find(y)

    for (ID, u, v, w) in G.E
        (Selected[u] && Selected[v]) ? SetUnion(u, v) : nothing
    end
    status = zeros(Int, G.n)

    while true
        clist = Array{Int, 1}()
        for i = 1 : G.n
            if Selected[i] && (F[i] == i)
                push!(clist, i)
            end
        end
        Bst = Array{Int, 1}()
        bstlen = G.n
        fill!(status, 0)
        for u in clist
            for i = 1 : G.n
                if Selected[i]
                    status[i] = (find(i) == u) ? 1 : 2
                end
            end
            tmp = SearchWay(g, status)
            if (tmp != nothing) && (size(tmp, 1) < bstlen)
                bstlen = size(tmp, 1)
                Bst = copy(tmp)
            end
        end
        if bstlen == G.n
            break
        end
        for i = 1 : bstlen
            Selected[Bst[i]] = true
            (i > 1) ? SetUnion(Bst[i], Bst[i - 1]) : nothing
        end
    end

    NodeC = Dict{Int, Float64}()
    NCmax = NC[argmax(NC)]
    for i = 1 : G.n
        if Selected[i]
            NodeC[G.V[i]] = NC[i] / NCmax
        end
    end

    return SubGraph(G, Selected, EC), NodeC, NCmax
end

function SearchAll(g, st, selected)
    n = size(g, 1)
    visited = zeros(Bool, n)
    visited[st] = true
    Q = zeros(Int, n)
    front = 1
    rear = 0

    pushQ(x) = begin
        rear += 1
        Q[rear] = x
    end
    popQ() = begin
        tmp = Q[front]
        front += 1
        return tmp
    end

    pushQ(st)

    prior = zeros(Int, n)

    while front <= rear
        u = popQ()
        for v in g[u]
            if !visited[v]
                visited[v] = true
                prior[v] = u
                pushQ(v)
            end
        end
    end

    needNodes = zeros(Bool, n)
    for v = 1 : n
        if selected[v]
            x = v
            while (x > 0)
                needNodes[x] = true
                x = prior[x]
            end
        end
    end

    return needNodes
end

function KeepAll(G, Selected, NC, EC)
    g = getAdjacentList(G)
    Needed = copy(Selected)
    for i = 1 : G.n
        if Selected[i]
            needNodes = SearchAll(g, i, Selected)
            foreach(j -> Needed[j] = Needed[j] || needNodes[j], 1 : G.n)
        end
    end

    NodeC = Dict{Int, Float64}()
    NCmax = NC[argmax(NC)]
    for i = 1 : G.n
        if Needed[i]
            NodeC[G.V[i]] = NC[i] / NCmax
        end
    end

    return SubGraph(G, Needed, EC), NodeC, NCmax
end

function Compress(G, selectedNodes, selectedEdges, selectMode, addMode)
    NC, EC, Selected = selectNodes(G, selectedNodes, selectedEdges, selectMode)
    return (addMode == "KeepOne") ? KeepOne(G, Selected, NC, EC) : KeepAll(G, Selected, NC, EC)
end

function CentralityInformation(G, NodeC, maxNodeC)
    nci = 0.0
    eci = 0.0
    foreach(i -> nci += (NodeC[G.V[i]]*maxNodeC), 1 : G.n)
    foreach(i -> eci += G.E[i][4], 1 : G.m)
    return nci, eci, (nci+eci)/2
end
