include("Graph.jl")
include("Planarity.jl")
include("VirtualEdge.jl")
include("PlanarDrawing.jl")

function DrawingAnalyze(G, Emb, X, Y, PATH; picSize = 20.0)
	ll = picSize / max(X[argmax(X)], Y[argmax(Y)])
	minAngle = 1e10

	vectorDot(p1, p2) = begin
		return p1[1]*p2[1] + p1[2]*p2[2]
	end
	vectorLen(p) = begin
		return sqrt(p[1]*p[1] + p[2]*p[2])
	end
	getNextIndex(x, y) = begin
		ret = (y == size(Emb[x], 1)) ? 1 : (y+1)
		while Emb[x][ret] > G.m
			ret = (ret == size(Emb[x], 1)) ? 1 : (ret+1)
		end
		return ret
	end
	getVector(u, eid) = begin
		if (X[u], Y[u]) == first(PATH[eid])
			return PATH[eid][2][1]-X[u], PATH[eid][2][2]-Y[u]
		end
		return PATH[eid][size(PATH[eid], 1)-1][1]-X[u], PATH[eid][size(PATH[eid], 1)-1][2]-Y[u]
	end
	getAngle(x, y) = begin
		if (y == getNextIndex(x, y))
			return 1e10
		end
		p1 = getVector(x, Emb[x][y])
		p2 = getVector(x, Emb[x][getNextIndex(x, y)])
		return acosd(vectorDot(p1, p2) / (vectorLen(p1) * vectorLen(p2)))
	end

	for i = 1 : G.n
		for j = 1 : size(Emb[i], 1)
			if Emb[i][j] <= G.m
				minAngle = min(minAngle, getAngle(i, j))
			end
		end
	end

	println("Minimum Angle : ", minAngle)

	getDistance(u, v) = begin
		return vectorLen((X[u]-X[v], Y[u]-Y[v]))
	end

	minP2Pdistance = 1e10
	for i = 1 : G.n, j = i+1 : G.n
		minP2Pdistance = min(minP2Pdistance, getDistance(i, j))
	end

	println("Minimum point-to-point distance : ", minP2Pdistance*ll)

	getiv(xa, ya, xb, yb, xx, yy) = begin
		lenab = vectorLen((xb-xa, yb-ya))
		ex = (xb - xa) / lenab
		ey = (yb - ya) / lenab
		vx = xx - xa
		vy = yy - ya
		pt = vx*ex + vy*ey
		xi = xa + ex*pt
		yi = ya + ey*pt
		return xi, yi
	end
	ivInAB(xi, yi, xa, ya, xb, yb) = begin
		if (xi < min(xa, xb)) || (xi > max(xa, xb))
			return false
		elseif (yi < min(ya, yb)) || (yi > max(ya, yb))
			return false
		end
		return true
	end
	getP2Ldistance(p, la, lb) = begin
		xi, yi = getiv(la[1], la[2], lb[1], lb[2], p[1], p[2])
		if !ivInAB(xi, yi, la[1], la[2], lb[1], lb[2])
			return 1e10
		end
		return vectorLen((p[1]-xi, p[2]-yi))
	end

	minP2Ldistance = 1e10

	updateP2L(x) = begin
		if (x < minP2Ldistance)
			minP2Ldistance = x
		end
	end

	for i = 1 : G.n, j = 1 : G.m
		if (G.E[j][2] != i) && (G.E[j][3] != i)
			for k = 1 : size(PATH[j], 1)-1
				updateP2L(getP2Ldistance((X[i], Y[i]), PATH[j][k], PATH[j][k+1]))
			end
		end
	end

	println("Minimum point-to-line distance : ", minP2Ldistance*ll)
end

function PlanarPolyDraw(G)
	# get planar subgraph
	edgeIndex = zeros(Int, G.m)
	foreach(i -> edgeIndex[i] = i, 1 : G.m)
	sort!(edgeIndex, by = x->-G.E[x][4])
	EdgeList = Array{Tuple{Int32, Int32, Int32, Float64}, 1}()
	preservedEdgeList = Array{Tuple{Int32, Int32}, 1}()
	deletedEdgeList = Array{Tuple{Int32, Int32}, 1}()
	deletedEdgeC = Array{Float64, 1}()
	preservedEdgeC = 0.0
	M2 = 0
	for i = 1 : G.m
		(ID, u, v, w) = G.E[edgeIndex[i]]
		M2 += 1
		push!(EdgeList, (M2, u, v, 1.0))
		if LRtest(Graph(G.n, size(EdgeList, 1), G.V, EdgeList)) == nothing
			M2 -= 1
			pop!(EdgeList)
			push!(deletedEdgeList, (u, v))
			push!(deletedEdgeC, w)
		else
			preservedEdgeC += w
			push!(preservedEdgeList, (u, v))
		end
	end

	# add some virtual edge to make the graph connected
	EdgeList = copy(AddToConnected(Graph(G.n, size(EdgeList, 1), G.V, EdgeList)))
	# add some virtual edge to make the graph biconnected
	Emb = copy(LRtest(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), needEmbedding = true))
	EdgeList = copy(AddToBiConnected(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), Emb))
	# add some virtual edge to get a maximal planar graph
	Emb = copy(LRtest(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), needEmbedding = true))
	EdgeList = copy(AddToMaximalPlanar(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), Emb))
	# calculate canonical ordering
	Emb = copy(LRtest(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), needEmbedding = true))
	Order = CanonicalOrdering(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), Emb)

	# using polyline drawing to draw a planar graph
	X, Y, PATH, L, H = PolylineDrawing(Graph(G.n, size(EdgeList, 1), G.V, EdgeList), M2, Order, Emb)

	#DrawingAnalyze(Graph(G.n, M2, G.V, EdgeList[1:M2]), Emb, X, Y, PATH)

	return X, Y, PATH, deletedEdgeList, L, H, preservedEdgeC, deletedEdgeC
end
