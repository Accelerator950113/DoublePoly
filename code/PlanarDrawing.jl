include("Graph.jl")

function CanonicalOrdering(G, Emb)
	processedNeighbor = zeros(Int, G.n)
	status = ones(Int, G.n)
	# status : 0 -> Processed
	#          1 -> Unprocessed
	#          2 -> One Neighbor Processed
	#          10 -> Ready To Be Processed
	Q = zeros(Int, 2*G.n)
	Order = Array{Int, 1}()
	front = 1
	rear = 0

	target(eid, u) = begin
		return (G.E[eid][2] != u) ? G.E[eid][2] : G.E[eid][3]
	end
	PushQ(x) = begin
		rear += 1
		Q[rear] = x
	end
	PopQ() = begin
		tmp = Q[front]
		front += 1
		return tmp
	end

	firstVertex = 1
	secondVertex = target(Emb[firstVertex][1], firstVertex)

	PushQ(firstVertex)
	status[firstVertex] = 10
	PushQ(secondVertex)
	status[secondVertex] = 10

	while front <= rear
		u = PopQ()
		if ((status[u] != 10) && (u != secondVertex))
			continue
		end
		for i = 1 : size(Emb[u], 1)
			priorEdge = (i == 1) ? Emb[u][size(Emb[u], 1)] : Emb[u][i-1]
			nextEdge = (i == size(Emb[u], 1)) ? Emb[u][1] : Emb[u][i+1]
			v = target(Emb[u][i], u)
			priorVertex = target(priorEdge, u)
			nextVertex = target(nextEdge, u)

			if status[v] == 1
				status[v] = 2
				processedNeighbor[v] = u
			elseif status[v] == 2
				x = processedNeighbor[v]
				if (((nextVertex == x) && !((firstVertex == u) && (secondVertex == x))) || ((priorVertex == x) && !((firstVertex == x) && (secondVertex == u))))
					status[v] = 10
				else
					status[v] = 10+1
				end
			elseif status[v] > 2
				processedBefore = (status[priorVertex] == 0)
				processedAfter = (status[nextVertex] == 0)
				if ((processedBefore == false) && (processedAfter == false))
					status[v] += 1
				end
				if ((processedBefore == true) && (processedAfter == true))
					status[v] -= 1
				end
			end
			if status[v] == 10
				PushQ(v)
			end
		end
		status[u] = 0
		push!(Order, u)
	end

	return Order
end

function PolylineDrawing(G, m, Order, Emb)
	edges = Array{Tuple{Int, Int}, 1}()
	A = zeros(Int32, G.n, G.n)
	newID = zeros(Int, G.n)
	foreach(i -> newID[Order[i]] = i, 1 : G.n)
	emb = Array{Array{Int, 1}, 1}(undef, G.n)
	for i = 1 : G.n
		emb[i] = []
		for j = 1 : size(Emb[Order[i]], 1)
			push!(emb[i], Emb[Order[i]][j])
		end
	end
	for (ID, u, v, w) in G.E
		nu = newID[u]
		nv = newID[v]
		A[nu, nv] = ID
		A[nv, nu] = ID
		(nu > nv) ? push!(edges, (nv, nu)) : push!(edges, (nu, nv))
	end
	InList = Array{Array{Int, 1}, 1}(undef, G.n)
	foreach(i -> InList[i] = [], 1 : G.n)
	OutList = Array{Array{Int, 1}, 1}(undef, G.n)
	foreach(i -> OutList[i] = [], 1 : G.n)

	target(u, eid) = (edges[eid][1] == u) ? edges[eid][2] : edges[eid][1]
	next(u, p) = (p == 1) ? size(emb[u], 1) : p-1
	prior(u, p) = (p == size(emb[u], 1)) ? 1 : p+1

	# Step 1 : calculate in and out list
	## calculate OutList[1]
	pos = 1
	while target(1, emb[1][pos]) != 2
		pos += 1
	end
	pos = next(1, pos)
	for i = 1 : size(emb[1], 1)
		if emb[1][pos] <= m
			push!(OutList[1], emb[1][pos])
		end
		pos = next(1, pos)
	end
	## calculate InList[2] and OutList[2]
	pos = 1
	while target(2, emb[2][pos]) != 1
		pos += 1
	end
	(emb[2][pos] <= m) ? push!(InList[2], emb[2][pos]) : nothing
	pos = next(2, pos)
	for i = 1 : size(emb[2], 1)-1
		if emb[2][pos] <= m
			push!(OutList[2], emb[2][pos])
		end
		pos = next(2, pos)
	end
	## calculate InList and OutList for node 3-n
	PRED = zeros(Int, G.n)
	SUCC = zeros(Int, G.n)
	SUCC[1] = 2
	PRED[2] = 1

	Val = zeros(Int, G.n)
	LL = zeros(Int, G.n)
	RR = zeros(Int, G.n)
	faceStart = 1
	sp = 2
	Val[1], Val[2] = 1, 2
	LL[1], LL[2] = 0, 1
	RR[1], RR[2] = 2, 0

	addBetween(lv, rv, x) = begin
		sp += 1
		LL[sp], Val[sp], RR[sp] = lv, x, rv
		RR[lv], LL[rv] = sp, sp
	end

	for k = 3 : G.n
		w = Array{Int, 1}()
		p = faceStart
		lv = -1
		rv = -1
		while p != 0
			if A[k, Val[p]] != 0
				push!(w, Val[p])
				(lv == -1) ? (lv = p) : nothing
				(lv != -1) ? (rv = p) : nothing
			end
			p = RR[p]
		end
		wn = size(w, 1)
		PRED[k] = w[1]
		SUCC[k] = w[wn]
		p = 1
		while target(k, emb[k][p]) != w[1]
			p += 1
		end
		q = next(k, p)
		while target(k, emb[k][q]) != w[wn]
			(emb[k][q] <= m) ? push!(OutList[k], emb[k][q]) : nothing
			q = next(k, q)
		end
		while true
			(emb[k][p] <= m) ? push!(InList[k], emb[k][p]) : nothing
			if target(k, emb[k][p]) == w[wn]
				break
			end
			p = prior(k, p)
		end
		addBetween(lv, rv, k)
	end

	indeg(x) = size(InList[x], 1)
	outdeg(x) = size(OutList[x], 1)

	# Step 2 : calculate outpoint and inpoint for each edge
	edgeinx = zeros(Int, m)
	edgeiny = zeros(Int, m)
	edgeoutx = zeros(Int, m)
	edgeouty = zeros(Int, m)
	outl = zeros(Int, G.n)
	outr = zeros(Int, G.n)
	inl = zeros(Int, G.n)
	inr = zeros(Int, G.n)

	isLinked(u, v) = ((u != 0) && (v != 0) && (A[u, v] <= m))

	for v = 1 : G.n
		## arrange outpoint
		outvPLUS = div(outdeg(v), 2)
		outvMINUS = div(outdeg(v)-1, 2)
		deltal = 1
		deltar = 1
		if indeg(v) >= 2
			outl[v], outr[v] = outvMINUS, outvPLUS
		elseif indeg(v) == 1
			if isLinked(v, SUCC[v])
				outl[v], outr[v] = outvPLUS, outvMINUS
				deltal = 0
			elseif isLinked(v, PRED[v])
				outl[v], outr[v] = outvMINUS, outvPLUS
				deltar = 0
			else
				outl[v], outr[v] = outvMINUS, outvPLUS
			end
		else
			outl[v], outr[v] = outvMINUS, outvPLUS
			deltal, deltar = 0, 0
		end
		tx, ty = -outl[v], deltal
		if outdeg(v) > 0
			for i = 1 : outl[v]
				edgeoutx[OutList[v][i]] = tx
				edgeouty[OutList[v][i]] = ty
				tx += 1
				ty += 1
			end
			edgeoutx[OutList[v][outl[v]+1]] = 0
			edgeouty[OutList[v][outl[v]+1]] = max(outl[v]+deltal-1, outr[v]+deltar-1, 0)
			tx, ty = 1, outr[v]+deltar-1
			for i = outl[v]+2 : outdeg(v)
				edgeoutx[OutList[v][i]] = tx
				edgeouty[OutList[v][i]] = ty
				tx += 1
				ty -= 1
			end
		end

		## arrange inpoint
		inl[v] = max(0, div(indeg(v)-3, 2))
		inr[v] = max(0, div(indeg(v)-2, 2))
		if indeg(v) >= 3
			edgeinx[InList[v][1]] = -inl[v]
			edgeiny[InList[v][1]] = 0
			tx, ty = -inl[v], -1
			for i = 2 : inl[v]+1
				edgeinx[InList[v][i]] = tx
				edgeiny[InList[v][i]] = ty
				tx += 1
				ty -= 1
			end
			edgeinx[InList[v][inl[v]+2]] = 0
			edgeiny[InList[v][inl[v]+2]] = -inr[v]
			tx, ty = 1, -inr[v]
			for i = inl[v]+3 : indeg(v)-1
				edgeinx[InList[v][i]] = tx
				edgeiny[InList[v][i]] = ty
				tx += 1
				ty += 1
			end
			edgeinx[InList[v][indeg(v)]] = inr[v]
			edgeiny[InList[v][indeg(v)]] = 0
		end
	end

	# Step 3 : arrange all nodes
	fill!(Val, 0)
	fill!(LL, 0)
	fill!(RR, 0)
	Father = zeros(Int, G.n)
	faceStart = 1
	sp = 2
	Val[1], Val[2] = 1, 2
	LL[1], LL[2] = 0, 1
	RR[1], RR[2] = 2, 0

	X = zeros(Int, G.n)
	Y = zeros(Int, G.n)
	X[1] = outl[1]
	X[2] = outr[1] + outl[2] + 1

	for k = 3 : G.n
		w = Array{Int, 1}()
		p = faceStart
		lv = -1
		rv = -1
		while p != 0
			if A[k, Val[p]] != 0
				push!(w, Val[p])
				(lv == -1) ? (lv = p) : nothing
				(lv != -1) ? (rv = p) : nothing
			end
			p = RR[p]
		end
		wn = size(w, 1)
		maxy = Y[w[1]]
		foreach(i -> maxy = max(maxy, Y[w[i]]), 2 : wn)
		Y[k] = inr[k] + maxy + 1
		sumx = 0
		for i = 2 : wn
			sumx += X[w[i]]
			X[w[i]] = sumx
		end
		## calculate nextLeftK and nextRightK
		nextLeftK, nextRightK = 0, 0
		pp = -1
		if outdeg(PRED[k]) > 0
			for i = 1 : outdeg(PRED[k])
				if target(PRED[k], OutList[PRED[k]][i]) == k
					pp = OutList[PRED[k]][i]
					break
				end
			end
			nextRightK = (pp == -1) ? (outr[PRED[k]]+1) : (edgeoutx[pp])
		end
		pp = -1
		if outdeg(SUCC[k]) > 0
			for i = 1 : outdeg(SUCC[k])
				if target(SUCC[k], OutList[SUCC[k]][i]) == k
					pp = OutList[SUCC[k]][i]
					break
				end
			end
			nextLeftK = (pp == -1) ? (-outl[SUCC[k]]-1) : (edgeoutx[pp])
		end
		dxl = nextRightK
		dxr = nextLeftK

		if indeg(k) >= 3
			dxt = edgeoutx[InList[k][inl[k]+2]]
			ct = target(k, InList[k][inl[k]+2])
			X[k] = max(X[ct] + dxt, dxl + outl[k])
			dta = X[k] - (X[ct] + dxt)
			X[SUCC[k]] = max(X[SUCC[k]] + dta - X[k], outr[k] - dxr)
			posT = 1
			while w[posT] != ct
				posT += 1
			end
			foreach(i -> X[w[i]] -= X[k], 2 : posT-1)
			foreach(i -> X[w[i]] += (dta - X[k]), posT : wn-1)
			foreach(i -> Father[w[i]] = k, 2 : wn-1)
		else
			X[k] = outl[k] + dxl
			(indeg(k) == 2 && isLinked(k, PRED[k]) == false) ? (X[k] = max(X[k], X[target(k, InList[k][1])] + edgeoutx[InList[k][1]])) : nothing
			if (indeg(k) == 2) && (!isLinked(k, SUCC[k])) && (X[k] > (X[target(k, InList[k][2])] + edgeoutx[InList[k][2]]))
				tmpdx = X[k] - (X[target(k, InList[k][2])] + edgeoutx[InList[k][2]])
				tmppt = 1
				while (w[tmppt] != target(k, InList[k][2]))
					tmppt += 1
				end
				foreach(i -> X[w[i]] += tmpdx, tmppt : wn)
			end
			X[SUCC[k]] = max(outr[k]-dxr, X[SUCC[k]]-X[k])
			foreach(i -> X[w[i]] -= X[k], 2 : wn-1)
			foreach(i -> Father[w[i]] = k, 2 : wn-1)
		end
		addBetween(lv, rv, k)
	end

	X2 = zeros(Int, G.n)
	finalw = Array{Int, 1}()
	p2 = faceStart
	while p2 != 0
		push!(finalw, Val[p2])
		p2 = RR[p2]
	end
	visited = zeros(Bool, G.n)
	sumx2 = 0
	for i = 1 : size(finalw, 1)
		sumx2 += X[finalw[i]]
		X2[finalw[i]] = sumx2
		visited[finalw[i]] = true
	end

	getX2(u) = begin
		if visited[u] == false
			X2[u] = X[u] + getX2(Father[u])
			visited[u] = true
		end
		return X2[u]
	end

	foreach(i -> X2[i] = getX2(i), 1 : G.n)

	# Step 4 : calculate path
	realX = zeros(Int, G.n)
	realY = zeros(Int, G.n)
	PATH = Array{Array{Tuple{Int, Int}, 1}, 1}(undef, m)
	deg = zeros(Int, G.n)
	foreach(i -> deg[i] = indeg(i) + outdeg(i), 1 : G.n)

	cross(x1, y1, x2, y2) = (x1*y2 - x2*y1)

	addToPATH(eid, x, y) = begin
		tn = size(PATH[eid], 1)
		if tn == 0
			push!(PATH[eid], (x, y))
		elseif (x != PATH[eid][tn][1]) || (y != PATH[eid][tn][2])
			if tn > 1
				if cross(PATH[eid][tn-1][1] - PATH[eid][tn][1], PATH[eid][tn-1][2] - PATH[eid][tn][2], x - PATH[eid][tn][1], y - PATH[eid][tn][2]) == 0
					PATH[eid][tn] = (x, y)
				else
					push!(PATH[eid], (x, y))
				end
			else
				push!(PATH[eid], (x, y))
			end
		end
	end

	for i = 1 : m
		(u, v) = edges[i]
		PATH[i] = []
		addToPATH(i, X2[u], Y[u])
		tx1, ty1 = X2[u]+edgeoutx[i], Y[u]+edgeouty[i]
		addToPATH(i, tx1, ty1)
		tx2, ty2 = X2[v]+edgeinx[i], Y[v]+edgeiny[i]
		((tx1 != tx2) && (ty1 != ty2)) ? addToPATH(i, tx1, ty2) : nothing
		addToPATH(i, tx2, ty2)
		addToPATH(i, X2[v], Y[v])
	end

	## dealing with one degree node
	hv = zeros(Bool, G.n)
	Snode = Array{Int, 1}()
	Sedge = Array{Int, 1}()
	Spath = Array{Tuple{Int, Int}, 1}()

	getEdgeID(u, pred) = begin
		for x in InList[u]
			if target(u, x) != pred
				return x, target(u, x)
			end
		end
		for x in OutList[u]
			if target(u, x) != pred
				return x, target(u, x)
			end
		end
		return nothing
	end

	addToSpath(u, eu) = begin
		if u < target(u, eu)
			for i = 1 : size(PATH[eu], 1)
				push!(Spath, PATH[eu][i])
			end
		else
			for i = size(PATH[eu], 1) : -1 : 1
				push!(Spath, PATH[eu][i])
			end
		end
	end

	getDir(x0, y0, x1, y1) = begin
		if (x0 == x1) && (y0 > y1)
			return 0, -1
		elseif (x0 == x1) && (y0 < y1)
			return 0, 1
		elseif (x0 > x1) && (y0 == y1)
			return -1, 0
		else
			return 1, 0
		end
	end

	DEC() = begin
		top = size(Spath, 1)
		while (top > 1) && (Spath[top] == Spath[top-1])
			pop!(Spath)
			top -= 1
		end
		if top > 1
			dx, dy = getDir(Spath[top][1], Spath[top][2], Spath[top-1][1], Spath[top-1][2])
			Spath[top] = (Spath[top][1]+dx, Spath[top][2]+dy)
		end
	end

	for i = 1 : G.n
		if (!hv[i]) && (deg[i] == 1)
			while size(Snode, 1) > 0
				pop!(Snode)
			end
			while size(Sedge, 1) > 0
				pop!(Sedge)
			end
			while size(Spath, 1) > 0
				pop!(Spath)
			end
			u = i
			eu, v = getEdgeID(u, -1)
			hv[v] = true
			push!(Snode, u)
			push!(Sedge, eu)
			addToSpath(u, eu)
			while deg[v] == 2
				predu = last(Snode)
				u = v
				eu, v = getEdgeID(u, predu)
				hv[v] = true
				push!(Snode, u)
				push!(Sedge, eu)
				addToSpath(u, eu)
			end
			tp = size(Spath, 1)
			# arrange the first node
			u = last(Snode)
			eu = last(Sedge)
			PATH[eu] = []
			push!(PATH[eu], last(Spath))
			if u < target(u, eu)
				if (edgeinx[eu] == 0) && (edgeiny[eu] == 0)
					DEC()
				else
					pop!(Spath)
					push!(Spath, (X2[target(u, eu)]+edgeinx[eu], Y[target(u, eu)]+edgeiny[eu]))
				end
			else
				if (edgeoutx[eu] == 0) && (edgeouty[eu] == 0)
					DEC()
				else
					pop!(Spath)
					push!(Spath, (X2[target(u, eu)]+edgeoutx[eu], Y[target(u, eu)]+edgeouty[eu]))
				end
			end
			push!(PATH[eu], last(Spath))
			X2[u] = last(Spath)[1]
			Y[u] = last(Spath)[2]
			pop!(Snode)
			pop!(Sedge)
			# arrange rest of the nodes
			while size(Snode, 1) > 0
				u = last(Snode)
				eu = last(Sedge)
				PATH[eu] = []
				push!(PATH[eu], last(Spath))
				DEC()
				push!(PATH[eu], last(Spath))
				X2[u] = last(Spath)[1]
				Y[u] = last(Spath)[2]
				pop!(Snode)
				pop!(Sedge)
			end
		end
	end

	## dealing with two degree node
	getTwoPaths(u) = begin
		ids = Vector{Int}()
		foreach(x -> push!(ids, x), InList[u])
		foreach(x -> push!(ids, x), OutList[u])
		for i in ids
			(PATH[i][1] != (X2[u], Y[u])) ? reverse!(PATH[i]) : nothing
		end
		return ids;
	end

	isBend(p1, p2, p3) = !(((p1[1] == p2[1]) && (p1[1] == p3[1])) || ((p1[2] == p2[2]) && (p1[2] == p3[2])))

	closePoint(u, v, eid) = begin
		tp = (u < v) ? (edgeoutx[eid], edgeouty[eid]) : (edgeinx[eid], edgeiny[eid])
		if (tp != (0, 0))
			return X2[u]+tp[1], Y[u]+tp[2]
		end
		dx, dy = getDir(X2[u], Y[u], X2[v], Y[v])
		return X2[u]+dx, Y[u]+dy
	end

	for i = 1 : G.n
		if (!hv[i]) && (deg[i] == 2)
			ids = getTwoPaths(i)
			if isBend(PATH[ids[1]][2], PATH[ids[1]][1], PATH[ids[2]][2])
				continue
			end
			if (size(PATH[ids[1]], 1) > 2)
				PATH[ids[1]] = PATH[ids[1]][2 : size(PATH[ids[1]], 1)]
				X2[i], Y[i] = PATH[ids[1]][1]
				PATH[ids[2]][1] = PATH[ids[1]][1]
			else
				X2[i], Y[i] = closePoint(target(i, ids[1]), i, ids[1])
				PATH[ids[1]][1] = (X2[i], Y[i])
				PATH[ids[2]][1] = (X2[i], Y[i])
			end
		end
	end

	# delete unused rows and columns
	Xs = Array{Int, 1}()
	Ys = Array{Int, 1}()
	NewX = Dict{Int, Int}()
	NewY = Dict{Int, Int}()
	for xx in X2
		push!(Xs, xx)
	end
	for yy in Y
		push!(Ys, yy)
	end
	for i = 1 : size(PATH, 1)
		for (xx, yy) in PATH[i]
			push!(Xs, xx)
			push!(Ys, yy)
		end
	end
	sort!(Xs)
	sort!(Ys)
	totalX = 1
	totalY = 1
	for xx in Xs
		if !haskey(NewX, xx)
			NewX[xx] = totalX
			totalX += 1
		end
	end
	for yy in Ys
		if !haskey(NewY, yy)
			NewY[yy] = totalY
			totalY += 1
		end
	end
	realPATH = Array{Array{Tuple{Int, Int}, 1}, 1}(undef, m)
	for i = 1 : size(PATH, 1)
		realPATH[i] = []
		for (xx, yy) in PATH[i]
			push!(realPATH[i], (NewX[xx], NewY[yy]))
		end
	end

	for i = 1 : G.n
		realX[Order[i]] = NewX[X2[i]]
		realY[Order[i]] = NewY[Y[i]]
	end
	return realX, realY, realPATH, totalX-1, totalY-1
end
