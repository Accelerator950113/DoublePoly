include("Graph.jl")

function AddLine(X, Y, PATH, edges, L, H)
    n = size(X, 1)
    m = size(edges, 1)
    leftX = ones(Int, L)
    rightX = ones(Int, L)
    upY = ones(Int, H)
    downY = ones(Int, H)

    deg = zeros(Int, n)
    for (u, v) in edges
        deg[u] += 1
        deg[v] += 1
    end
    for i = 1 : n
        x, y = X[i], Y[i]
        leftX[x] = max(leftX[x], div(deg[i]+3, 4)+1)
        rightX[x] = max(rightX[x], div(deg[i]+3, 4)+1)
        upY[y] = max(upY[y], div(deg[i]+3, 4)+1)
        downY[y] = max(downY[y], div(deg[i]+3, 4)+1)
    end

    dicX = Dict{Int, Bool}()
	dicY = Dict{Int, Bool}()
    nX = zeros(Int, L)
    nY = zeros(Int, H)
    nL, nH = 0, 0
    for i = 1 : L
        nL += (i == 1) ? (leftX[i]+1) : (leftX[i]+rightX[i-1]+1)
        nX[i] = nL
		dicX[nL] = true
    end
    nL += rightX[L]
    for i = 1 : H
        nH += (i == 1) ? (downY[i]+1) : (downY[i]+upY[i-1]+1)
        nY[i] = nH
		dicY[nL] = true
    end
    nH += upY[H]

    MAP = zeros(Bool, nL, nH)
    used = zeros(Bool, nL, nH)
    for i = 1 : n
        MAP[nX[X[i]], nY[Y[i]]] = true
        rr = div(deg[i]+3, 4)
        if rr > 0
            xx, yy = nX[X[i]], nY[Y[i]]
            di, dj = rr, 0
            while dj < rr
                MAP[xx+di, yy+dj] = true
                MAP[xx-dj, yy+di] = true
                MAP[xx+dj, yy-di] = true
                MAP[xx-di, yy-dj] = true
                di -= 1
                dj += 1
            end
        end
    end
    dx = [0, 0, 1, -1]
    dy = [1, -1, 0, 0]
    isInQ = zeros(Int32, nL, nH)

	cross(x1, y1, x2, y2) = (x1*y2 - x2*y1)
	isBend(x1, y1, x2, y2, x3, y3) = (cross(x1-x2, y1-y2, x3-x2, y3-y2) != 0)

    # find a shortest path from node u to node v
    SearchWay(u, v) = begin
        target = Dict{Tuple{Int, Int}, Bool}()
        Q = Array{Tuple{Int, Int}, 1}()
        fill!(isInQ, 0)
        prior = Array{Int, 1}()
		bend = Array{Int, 1}()
        front, rear = 1, 0
        pushQ(x, y) = begin
            rear += 1
            push!(Q, (x, y))
            isInQ[x, y] = rear
            push!(prior, 0)
			push!(bend, 0)
        end
        popQ() = begin
            tp = Q[front]
            front += 1
            return tp
        end

        ## push start points into the queue
        rr = div(deg[u]+3, 4)
        xx, yy = nX[X[u]], nY[Y[u]]
        di, dj = rr, 0
        while dj < rr
            (used[xx+di, yy+dj] == false) ? pushQ(xx+di, yy+dj) : nothing
            (used[xx-dj, yy+di] == false) ? pushQ(xx-dj, yy+di) : nothing
            (used[xx+dj, yy-di] == false) ? pushQ(xx+dj, yy-di) : nothing
            (used[xx-di, yy-dj] == false) ? pushQ(xx-di, yy-dj) : nothing
            di -= 1
            dj += 1
        end

        ## mark target points
        rr = div(deg[v]+3, 4)
        xx, yy = nX[X[v]], nY[Y[v]]
        di, dj = rr, 0
        while dj < rr
            (used[xx+di, yy+dj] == false) ? (target[(xx+di, yy+dj)] = true) : nothing
            (used[xx-dj, yy+di] == false) ? (target[(xx-dj, yy+di)] = true) : nothing
            (used[xx+dj, yy-di] == false) ? (target[(xx+dj, yy-di)] = true) : nothing
            (used[xx-di, yy-dj] == false) ? (target[(xx-di, yy-dj)] = true) : nothing
            di -= 1
            dj += 1
        end

        checkRange(x, y) = ((x > 0) && (y > 0) && (x <= nL) && (y <= nH))

        ans = Array{Tuple{Int, Int}, 1}()
		minBend = 1048576
        ## find the shortest path
        while front <= rear
            cntx, cnty = popQ()
			prex, prey = (prior[front-1] == 0) ? (nX[X[u]], nY[Y[u]]) : Q[prior[front-1]]
            for i = 1 : 4
                tx, ty = cntx+dx[i], cnty+dy[i]
				if (checkRange(tx, ty) == false) || (used[tx, ty] == true)
					continue
				end
				cntbend = isBend(prex, prey, cntx, cnty, tx, ty) ? (bend[front-1]+1) : bend[front-1]
				#cntbend += (!haskey(dicX, cntx)) + (!haskey(dicY, cnty))
                if haskey(target, (tx, ty)) && (cntbend+isBend(cntx, cnty, tx, ty, nX[X[v]], nY[Y[v]]) < minBend)
					minBend = cntbend + isBend(cntx, cnty, tx, ty, nX[X[v]], nY[Y[v]]) 
					tans = Array{Tuple{Int, Int}, 1}()
                    push!(tans, (nX[X[v]], nY[Y[v]]))
                    push!(tans, (tx, ty))
                    pt = front-1
                    while pt != 0
                        push!(tans, Q[pt])
                        pt = prior[pt]
                    end
                    push!(tans, (nX[X[u]], nY[Y[u]]))
                    ans = copy(tans)
                end
                if (MAP[tx, ty] == false) && (isInQ[tx, ty] == 0)
                    pushQ(tx, ty)
                    prior[rear] = front-1
					bend[rear] = cntbend
                elseif (MAP[tx, ty] == false) && (isInQ[tx, ty] >= front)
					if bend[isInQ[tx, ty]] > cntbend
						bend[isInQ[tx, ty]] = cntbend
						prior[isInQ[tx, ty]] = front-1
					end
				end
            end
        end
        return (minBend == 1048576) ? nothing : ans
    end

    PATH2 = Array{Array{Tuple{Int, Int}, 1}, 1}(undef, m)

	addToPATH2(eid, x, y) = begin
		tn = size(PATH2[eid], 1)
		if tn == 0
			push!(PATH2[eid], (x, y))
		elseif (x != PATH2[eid][tn][1]) || (y != PATH2[eid][tn][2])
			if tn > 1
				if cross(PATH2[eid][tn-1][1] - PATH2[eid][tn][1], PATH2[eid][tn-1][2] - PATH2[eid][tn][2], x - PATH2[eid][tn][1], y - PATH2[eid][tn][2]) == 0
					PATH2[eid][tn] = (x, y)
				else
					push!(PATH2[eid], (x, y))
				end
			else
				push!(PATH2[eid], (x, y))
			end
		end
	end

    dropedEdges = 0
    for i = 1 : m
        tmp = SearchWay(edges[i][1], edges[i][2])
        PATH2[i] = []
		if tmp == nothing
			dropedEdges += 1
			continue
		end
		for (xx, yy) in tmp
			addToPATH2(i, xx, yy)
			used[xx, yy] = true
		end
    end

	Xs = Array{Int, 1}()
	Ys = Array{Int, 1}()
	NewX = Dict{Int, Int}()
	NewY = Dict{Int, Int}()
	totalX, totalY = 0, 0
	for i = 1 : n
		push!(Xs, nX[X[i]])
		push!(Ys, nY[Y[i]])
	end
	for i = 1 : size(PATH, 1)
		for (xx, yy) in PATH[i]
			push!(Xs, nX[xx])
			push!(Ys, nY[yy])
		end
	end
	for i = 1 : size(PATH2, 1)
		for (xx, yy) in PATH2[i]
			push!(Xs, xx)
			push!(Ys, yy)
		end
	end
	sort!(Xs)
	sort!(Ys)
	for xx in Xs
		if !haskey(NewX, xx)
			totalX += 1
			NewX[xx] = totalX
		end
	end
	for yy in Ys
		if !haskey(NewY, yy)
			totalY += 1
			NewY[yy] = totalY
		end
	end

	realX = zeros(Int, n)
	realY = zeros(Int, n)
	realPATH = Array{Array{Tuple{Int, Int}, 1}, 1}(undef, size(PATH, 1))
	realPATH2 = Array{Array{Tuple{Int, Int}, 1}, 1}(undef, m)
	for i = 1 : n
		realX[i] = NewX[nX[X[i]]]
		realY[i] = NewY[nY[Y[i]]]
	end
	for i = 1 : size(PATH, 1)
		realPATH[i] = []
		for (xx, yy) in PATH[i]
			push!(realPATH[i], (NewX[nX[xx]], NewY[nY[yy]]))
		end
	end
	for i = 1 : size(PATH2, 1)
		realPATH2[i] = []
		for (xx, yy) in PATH2[i]
			push!(realPATH2[i], (NewX[xx], NewY[yy]))
		end
	end

	# reuse rows and columns
	statusX = zeros(Int8, totalX)
	statusY = zeros(Int8, totalY)
	for i = 1 : n
		statusX[realX[i]] |= 3
		statusY[realY[i]] |= 3
	end
	for i = 1 : size(realPATH, 1)
		for (xx, yy) in realPATH[i]
			statusX[xx] |= 1
			statusY[yy] |= 1
		end
	end
	for i = 1 : size(realPATH2, 1)
		for (xx, yy) in realPATH2[i]
			statusX[xx] |= 2
			statusY[yy] |= 2
		end
	end

	X2 = zeros(Int, totalX)
	Y2 = zeros(Int, totalY)
	ansX1 = 0
	ansX2 = 0
	ansY1 = 0
	ansY2 = 0

	for i = 1 : totalX
		if statusX[i] == 3
			tmp = max(ansX1, ansX2)
			ansX1 = tmp+1
			ansX2 = tmp+1
			X2[i] = tmp+1
		elseif statusX[i] == 2
			ansX2 += 1
			X2[i] = ansX2
		else
			ansX1 += 1
			X2[i] = ansX1
		end
	end

	for i = 1 : totalY
		if statusY[i] == 3
			tmp = max(ansY1, ansY2)
			ansY1 = tmp+1
			ansY2 = tmp+1
			Y2[i] = tmp+1
		elseif statusY[i] == 2
			ansY2 += 1
			Y2[i] = ansY2
		else
			ansY1 += 1
			Y2[i] = ansY1
		end
	end
	map!(x -> X2[x], realX, realX)
	map!(y -> Y2[y], realY, realY)
	foreach(i -> map!(p -> (X2[p[1]], Y2[p[2]]), realPATH[i], realPATH[i]), 1 : size(realPATH, 1))
	foreach(i -> map!(p -> (X2[p[1]], Y2[p[2]]), realPATH2[i], realPATH2[i]), 1 : size(realPATH2, 1))
	totalX = max(ansX1, ansX2)
	totalY = max(ansY1, ansY2)

	return realX, realY, realPATH, realPATH2, totalX, totalY, dropedEdges
end

function CalculatePreservedCentralityInformation(G, NodeC, maxNodeC, preservedEdgeC, deletedEdgeC, PATH2)
	nodeInfo = 0.0
	edgeInfo = preservedEdgeC
	for i = 1 : G.n
		nodeInfo += (NodeC[G.V[i]]*maxNodeC)
	end
	for i = 1 : size(PATH2, 1)
		if size(PATH2[i], 1) > 0
			edgeInfo += deletedEdgeC[i]
		end
	end

	return nodeInfo, edgeInfo, (nodeInfo+edgeInfo)/2
end
