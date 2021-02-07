struct Graph
    n :: Int32 # |V|
    m :: Int32 # |E|
    V :: Array{Int32, 1} # V[i] = Real Index of node i
    E :: Array{Tuple{Int32, Int32, Int32, Float64}, 1} # (ID, u, v, w) in Edge Set
end

function readGraph(fileName, graphType) # Read graph from a file | graphType = [weighted, unweighted]
    # Initialized
    n = 0
    origin = Dict{Int32, Int32}()
    label = Dict{Int32, Int32}()
    edge = Set{Tuple{Int32, Int32, Float64}}()

    getid(x :: Int32) = haskey(label, x) ? label[x] : label[x] = n += 1

    open(fileName) do f1
        for line in eachline(f1)
            # Read origin data from file
            buf = split(line)
            u = parse(Int32, buf[1])
            v = parse(Int32, buf[2])
            if graphType == "weighted"
                w = parse(Float64, buf[3])
            else
                w = 1.0
            end
            if u == v
                continue
            end
            # Label the node
            u1 = getid(u)
            v1 = getid(v)
            origin[u1] = u
            origin[v1] = v
            # Store the edge
            if u1 > v1
                u1, v1 = v1, u1
            end
            push!(edge, (u1, v1, w))
        end
    end

    # Store data into the struct Graph
    m = length(edge)
    V = Array{Int32, 1}(undef, n)
    E = Array{Tuple{Int32, Int32, Int32, Float64}, 1}(undef, m)

    for i = 1 : n
        V[i] = origin[i]
    end

    ID = 0
    for (u, v, w) in edge
        ID = ID + 1
        E[ID] = (ID, u, v, w)
    end

    return Graph(n, m, V, E)
end

function getConnectedComponents(G)
    F = zeros(Int, G.n)
    foreach(i -> F[i] = i, 1 : G.n)

    find(x) = begin
        if F[x] != x
            F[x] = find(F[x])
        end
        return F[x]
    end

    for (ID, u, v, w) in G.E
        p = find(u)
        q = find(v)
        (p != q) ? F[p] = q : nothing
    end

    CC = zeros(Int, G.n)
    foreach(i -> CC[i] = find(i), 1 : G.n)

    return CC
end

function getBiconnectedComponents(G)
    dfn = zeros(Int, G.n)
    low = zeros(Int, G.n)
    cnt = 0
    cutPoint = zeros(Bool, G.n)
    bccid = zeros(Int, G.n)
    edgeC = zeros(Int, G.m)
    S = zeros(Int, 2*G.m)
    top = 0

    nbcc = 0
    bcc = Array{Array{Int, 1}, 1}(undef, G.n)
    g = Array{Array{Int, 1}, 1}(undef, G.n)
    gID = Array{Array{Int, 1}, 1}(undef, G.n)
    for i = 1 : G.n
        g[i] = []
        gID[i] = []
        bcc[i] = []
    end
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(gID[u], ID)
        push!(g[v], u)
        push!(gID[v], ID)
    end

    DFS(u, fa) = begin
        cnt += 1
        dfn[u] = cnt
        low[u] = cnt
        son = 0
        for i = 1 : size(g[u], 1)
            v = g[u][i]
            if dfn[v] == 0
                top += 1
                S[top] = gID[u][i]
                son += 1
                DFS(v, u)
                low[u] = min(low[u], low[v])
                if dfn[u] <= low[v]
                    cutPoint[u] = true
                    nbcc += 1
                    while true
                        x = S[top]
                        top -= 1
                        if bccid[G.E[x][2]] != nbcc
                            bccid[G.E[x][2]] = nbcc
                            push!(bcc[nbcc], G.E[x][2])
                        end
                        if bccid[G.E[x][3]] != nbcc
                            bccid[G.E[x][3]] = nbcc
                            push!(bcc[nbcc], G.E[x][3])
                        end
                        if x == gID[u][i]
                            break
                        end
                    end
                end
            elseif ((v != fa) && (dfn[v] < dfn[u]))
                top += 1
                S[top] = gID[u][i]
                low[u] = min(low[u], dfn[v])
            end
        end
        if (fa == 0) && (son <= 1)
            cutPoint[u] = false
        end
    end

    for i = 1 : G.n
        if dfn[i] == 0
            DFS(i, 0)
        end
    end

    timeStamp = 0
    fill!(bccid, 0)
    for i = 1 : nbcc
        timeStamp += 1
        for j = 1 : size(bcc[i], 1)
            bccid[bcc[i][j]] = timeStamp
        end
        for u in bcc[i]
            for j = 1 : size(g[u], 1)
                v = g[u][j]
                ei = gID[u][j]
                (bccid[v] == timeStamp) ? edgeC[ei] = timeStamp : nothing
            end
        end
    end
    for i = 1 : G.m
        if edgeC[i] == 0
            timeStamp += 1
            edgeC[i] = timeStamp
        end
    end

    cutList = Array{Int, 1}()
    for i = 1 : G.n
        if cutPoint[i] == true
            push!(cutList, i)
        end
    end

    return cutList, edgeC
end

function IsBiconnected(G)
    cutList, edgeC = getBiconnectedComponents(G)
    return size(cutList, 1) == 0
end

function getAdjacentList(G)
    g = Array{Array{Int, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    return g
end

function SubGraph(G, Selected, EC)
    n = 0
    label = Dict{Int32, Int32}()
    nodeList = Array{Int32, 1}()

    for i = 1 : G.n
        if Selected[i]
            n += 1
            push!(nodeList, i)
            label[i] = n
        end
    end

    m = 0
    edges = Array{Tuple{Int32, Int32, Int32, Float64}, 1}()
    for (ID, u, v, w) in G.E
        if haskey(label, u) && haskey(label, v)
            m += 1
            push!(edges, (m, label[u], label[v], EC[ID]))
        end
    end

    V = zeros(Int, n)
    foreach(i -> V[i] = G.V[nodeList[i]], 1 : n)

    return Graph(n, m, V, edges)
end
