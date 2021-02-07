include("Graph.jl")

function LRtest(G; needEmbedding = false)
    if (G.n > 2) && (G.m > 3*G.n-6)
        return nothing
    end
    INF = 131072

    g = Array{Array{Int, 1}, 1}(undef, G.n)
    gID = Array{Array{Int, 1}, 1}(undef, G.n)
    for i = 1 : G.n
        g[i] = []
        gID[i] = []
    end
    for (ID, u, v, w) in G.E
        push!(g[u], v)
        push!(gID[u], ID)
        push!(g[v], u)
        push!(gID[v], ID)
    end

    Roots = Array{Int, 1}()
    height = zeros(Int, G.n)
    lowpt = zeros(Int, G.m)
    lowpt2 = zeros(Int, G.m)
    nestingDepth = zeros(Int, G.m)
    fill!(height, INF)
    parentEdge = zeros(Int, G.n)
    ori = zeros(Int, G.m)

    ref = zeros(Int, G.m)
    side = ones(Int, G.m)
    S = zeros(Int, 2*G.m, 4)
    top = 0
    stackBottom = zeros(Int, G.m)
    lowptEdge = zeros(Int, G.m)

    # orientation
    DFS1(v) = begin
        e = parentEdge[v]
        for i = 1 : size(g[v], 1)
            w = g[v][i]
            wid = gID[v][i]
            if ori[wid] == 0
                ori[wid] = (v == G.E[wid][2]) ? 1 : -1
                lowpt[wid] = height[v]
                lowpt2[wid] = height[v]
                if height[w] == INF
                    parentEdge[w] = wid
                    height[w] = height[v] + 1
                    DFS1(w)
                else
                    lowpt[wid] = height[w]
                end
                # determine nesting depth
                nestingDepth[wid] = 2*lowpt[wid]
                if lowpt2[wid] < height[v]
                    nestingDepth[wid] += 1
                end
                # update lowpoints of parent edge e
                if e != 0
                    if lowpt[wid] < lowpt[e]
                        lowpt2[e] = min(lowpt[e], lowpt2[wid])
                        lowpt[e] = lowpt[wid]
                    elseif lowpt[wid] > lowpt[e]
                        lowpt2[e] = min(lowpt2[e], lowpt[wid])
                    else
                        lowpt2[e] = min(lowpt2[e], lowpt2[wid])
                    end
                end
            end
        end
    end

    for i = 1 : G.n
        if height[i] == INF
            push!(Roots, i)
            height[i] = 0
            DFS1(i)
        end
    end

    #println("Lowpt :")
    #foreach(i -> print(lowpt[i], " "), 1 : G.m)
    #println()
    #println("Lowpt2 :")
    #foreach(i -> print(lowpt2[i], " "), 1 : G.m)
    #println()
    #println("Height :")
    #foreach(i -> print(height[i], " "), 1 : G.n)
    #println()
    #for i = 1 : G.m
    #    if ori[i] == 1
    #        println(G.E[i][2], " -> ", G.E[i][3])
    #    elseif ori[i] == -1
    #        println(G.E[i][3], " -> ", G.E[i][2])
    #    end
    #end

    # testing

    ep = Array{Array{Int, 1}, 1}(undef, G.n)
    eptar = Array{Array{Int, 1}, 1}(undef, G.n)
    foreach(i -> ep[i] = [], 1 : G.n)
    foreach(i -> eptar[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        if ori[ID] == 1
            push!(ep[u], ID)
        else
            push!(ep[v], ID)
        end
    end
    for i = 1 : G.n
        sort!(ep[i], by = x->nestingDepth[x])
        for j = 1 : size(ep[i], 1)
            if ori[ep[i][j]] == 1
                push!(eptar[i], G.E[ep[i][j]][3])
            else
                push!(eptar[i], G.E[ep[i][j]][2])
            end
        end
    end

    #for i = 1 : G.n
    #    for j = 1 : size(ep[i], 1)
    #        print(eptar[i][j], "[", ep[i][j], "](", nestingDepth[ep[i][j]], ") ")
    #    end
    #    println()
    #end

    conflicting(I, b) = begin
        return ((I != (0, 0)) && (lowpt[I[2]] > lowpt[b]))
    end
    lowest(P) = begin
        if (P[1] == 0) && (P[2] == 0)
            return lowpt[P[3]]
        elseif (P[3] == 0) && (P[4] == 0)
            return lowpt[P[1]]
        end
        return min(lowpt[P[1]], lowpt[P[3]])
    end
    target(x) = begin
        if ori[x] == 1
            return G.E[x][3]
        else
            return G.E[x][2]
        end
    end
    swapLR(Q) = begin
        return [Q[3], Q[4], Q[1], Q[2]]
    end

    DFS2(v) = begin
        e = parentEdge[v]
        for i = 1 : size(ep[v], 1)
            eiid = ep[v][i]
            eitar = eptar[v][i]
            stackBottom[eiid] = top
            if eiid == parentEdge[eitar]
                if !DFS2(eitar)
                    return false
                end
            else
                lowptEdge[eiid] = eiid
                top += 1
                S[top, :] = [0, 0, eiid, eiid]
            end
            # integrate new return edges
            if lowpt[eiid] < height[v]
                if i == 1
                    lowptEdge[e] = lowptEdge[eiid]
                else
                    # add constraints of ei
                    P = zeros(Int, 4)
                    # merge return edges of ei into P.R
                    while true
                        Q = copy(S[top, :])
                        top -= 1
                        ((Q[1] != 0) || (Q[2] != 0)) ? (Q = copy(swapLR(Q))) : nothing
                        if (Q[1] != 0) || (Q[2] != 0)
                            return false
                        else
                            if lowpt[Q[3]] > lowpt[e]
                                if (P[3] == 0) && (P[4] == 0)
                                    P[4] = Q[4]
                                else
                                    ref[P[3]] = Q[4]
                                end
                                P[3] = Q[3]
                            else
                                ref[Q[3]] = lowptEdge[e]
                            end
                        end
                        if (top == stackBottom[eiid])
                            break
                        end
                    end
                    # merge confliting return edges of e1, ..., ei-1 into P.L
                    #println(eiid, " ", S[top, 1], " ", S[top, 2], " ", S[top, 3], " ", S[top, 4])
                    while ((conflicting((S[top, 1], S[top, 2]), eiid)) || (conflicting((S[top, 3], S[top, 4]), eiid)))
                        Q = copy(S[top, :])
                        top -= 1
                        conflicting((Q[3], Q[4]), eiid) ? (Q = copy(swapLR(Q))) : nothing
                        if conflicting((Q[3], Q[4]), eiid)
                            return false
                        else
                            ref[P[3]] = Q[4]
                            (Q[3] != 0) ? (P[3] = Q[3]) : nothing
                        end
                        if (P[1] == 0) && (P[2] == 0)
                            P[2] = Q[2]
                        else
                            ref[P[1]] = Q[2]
                        end
                        P[1] = Q[1]
                    end
                    if sum(P) != 0
                        top += 1
                        S[top, :] = copy(P)
                    end
                end
            end
        end
        # remove back edges returning to parent
        if e != 0
            u = (G.E[e][2] != v) ? (G.E[e][2]) : (G.E[e][3])
            # trim back edges ending at parent u
            ## drop entire conflict pairs
            while (top > 0) && (lowest(S[top, :]) == height[u])
                P = copy(S[top, :])
                top -= 1
                (P[1] != 0) ? (side[P[1]] = -1) : nothing
            end
            if top > 0
                P = copy(S[top, :])
                top -= 1
                # trim left interval
                while (P[2] != 0) && (target(P[2]) == u)
                    P[2] = ref[P[2]]
                end
                if (P[2] == 0) && (P[1] != 0)
                    ref[P[1]] = P[3]
                    side[P[1]] = -1
                    P[1] = 0
                end
                # trim right interval
                while (P[4] != 0) && (target(P[4]) == u)
                    P[4] = ref[P[4]]
                end
                if (P[4] == 0) && (P[3] != 0)
                    ref[P[3]] = P[1]
                    side[P[3]] = -1
                    P[3] = 0
                end
                top += 1
                S[top, :] = copy(P)
            end
            # side of e is side of a highest return edge
            if lowpt[e] < height[u]
                hl = S[top, 2]
                hr = S[top, 4]
                if (hl != 0) && ((hr == 0) || (lowpt[hl] > lowpt[hr]))
                    ref[e] = hl
                else
                    ref[e] = hr
                end
            end
        end
        return true
    end

    for s in Roots
        if !DFS2(s)
            return nothing
        end
    end

    #foreach(i -> print(side[i], " "), 1 : G.m)
    #println()

    Emb = Array{Array{Int, 1}, 1}(undef, G.n)
    foreach(i -> Emb[i] = [], 1 : G.n)

    if (needEmbedding == false)
        return Emb
    end

    # get embedding

    sign(eid) = begin
        if ref[eid] != 0
            side[eid] = side[eid] * sign(ref[eid])
            ref[eid] = 0
        end
        return side[eid]
    end

    foreach(i -> nestingDepth[i] = sign(i) * nestingDepth[i], 1 : G.m)

    foreach(i -> ep[i] = [], 1 : G.n)
    foreach(i -> eptar[i] = [], 1 : G.n)
    for (ID, u, v, w) in G.E
        if ori[ID] == 1
            push!(ep[u], ID)
        else
            push!(ep[v], ID)
        end
    end
    for i = 1 : G.n
        sort!(ep[i], by = x->nestingDepth[x])
        for j = 1 : size(ep[i], 1)
            if ori[ep[i][j]] == 1
                push!(eptar[i], G.E[ep[i][j]][3])
            else
                push!(eptar[i], G.E[ep[i][j]][2])
            end
        end
    end

    sp = 0
    LL = zeros(Int, 3*G.m)
    RR = zeros(Int, 3*G.m)
    Val = zeros(Int, 3*G.m)
    Fst = zeros(Int, G.n)
    rightRef = zeros(Int, G.n)
    leftRef = zeros(Int, G.n)
    lkid = Array{Array{Int, 1}, 1}(undef, G.n)
    foreach(i -> lkid[i] = [], 1 : G.n)

    AddToFirst(v, x) = begin
        sp += 1
        RR[sp] = Fst[v]
        (Fst[v] != 0) ? (LL[Fst[v]] = sp) : nothing
        Val[sp] = x
        Fst[v] = sp
    end

    AddBefore(v, id, x) = begin
        sp += 1
        RR[sp] = id
        LL[sp] = LL[id]
        (LL[id] != 0) ? (RR[LL[id]] = sp) : nothing
        (Fst[v] == id) ? (Fst[v] = sp) : nothing
        LL[id] = sp
        Val[sp] = x
    end

    AddAfter(id, x) = begin
        sp += 1
        RR[sp] = RR[id]
        LL[sp] = id
        (RR[id] != 0) ? (LL[RR[id]] = sp) : nothing
        RR[id] = sp
        Val[sp] = x
    end

    for i = 1 : G.n
        if size(ep[i], 1) > 0
            AddToFirst(i, ep[i][1])
            push!(lkid[i], sp)
            for j = 2 : size(ep[i], 1)
                AddAfter(sp, ep[i][j])
                push!(lkid[i], sp)
            end
        end
    end

    DFS3(v) = begin
        for i = 1 : size(ep[v], 1)
            eiid = ep[v][i]
            w = eptar[v][i]
            eilkid = lkid[v][i]
            if eiid == parentEdge[w]
                AddToFirst(w, eiid)
                leftRef[v] = eilkid
                rightRef[v] = eilkid
                DFS3(w)
            else
                if side[eiid] == 1
                    AddAfter(rightRef[w], eiid)
                else
                    AddBefore(w, leftRef[w], eiid)
                    leftRef[w] = sp
                end
            end
        end
    end

    for s in Roots
        DFS3(s)
    end

    for i = 1 : G.n
        x = Fst[i]
        while x != 0
            push!(Emb[i], Val[x])
            x = RR[x]
        end
    end

    return Emb
end
