include("Graph.jl")

function AddToConnected(G)
	CC = getConnectedComponents(G)
	EdgeList = copy(G.E)
	nodes = zeros(Int, G.n)
	foreach(i -> nodes[i] = i, 1 : G.n)
	sort!(nodes, by = x->CC[x])
	M3 = G.m
	for i = 2 : G.n
		if CC[nodes[i]] != CC[nodes[i - 1]]
			M3 += 1
			push!(EdgeList, (M3, nodes[i], nodes[i - 1], 1.0))
		end
	end

	return EdgeList
end

function AddToBiConnected(G, Emb)
	cutList, edgeC = getBiconnectedComponents(G)
	EdgeList = copy(G.E)
	m = G.m

	target(u, eiid) = begin
		return (G.E[eiid][2] == u) ? (G.E[eiid][3]) : (G.E[eiid][2])
	end

	for u in cutList
		for i = 2 : size(Emb[u], 1)
			x = Emb[u][i-1]
			y = Emb[u][i]
			if edgeC[x] != edgeC[y]
				m += 1
				push!(EdgeList, (m, target(u, x), target(u, y), 1.0))
			end
		end
	end

	return EdgeList
end

function AddToMaximalPlanar(G, Emb)
	m = G.m
	nextEdge = zeros(Int, m, 2)
	visited = zeros(Bool, m, 2)
	EdgeList = copy(G.E)
	g = Array{Array{Int, 1}, 1}(undef, G.n)
	foreach(i -> g[i] = [], 1 : G.n)
	onFace = zeros(Int, G.n)
	timeStamp = 0

	for (ID, u, v, w) in G.E
		push!(g[u], v)
		push!(g[v], u)
	end

	getPos(eid, u) = begin
		return (G.E[eid][2] == u) ? 1 : 2
	end
	target(eid, u) = begin
		return (G.E[eid][2] != u) ? G.E[eid][2] : G.E[eid][3]
	end
	getNext(eid, u) = begin
		return (G.E[eid][2] == u) ? nextEdge[eid, 1] : nextEdge[eid, 2]
	end
	setVisited(eid, u) = begin
		if G.E[eid][2] == u
			visited[eid, 1] = true
		else
			visited[eid, 2] = true
		end
	end
	addEdge(u, v) = begin
		m += 1
		push!(EdgeList, (m, u, v, 1.0))
		push!(g[u], v)
		push!(g[v], u)
	end

	for u = 1 : G.n
		for i = 1 : size(Emb[u], 1)
			nextEdge[Emb[u][i], getPos(Emb[u][i], u)] = (i == size(Emb[u], 1)) ? Emb[u][1] : Emb[u][i+1]
		end
	end

	getFace(eid, u) = begin
		visited[eid, getPos(eid, u)] = true
		face = Array{Int, 1}()
		push!(face, u)
		cnte = eid
		cntu = u
		while true
			v = target(cnte, cntu)
			if v == face[1]
				break
			end
			push!(face, v)
			cnte = getNext(cnte, v)
			cntu = v
			visited[cnte, getPos(cnte, cntu)] = true
		end
		return face
	end

	m0 = m
	for i = 1 : m0
		for j = 1 : 2
			if visited[i, j] == false
				face = copy(getFace(i, G.E[i][j+1]))
				if size(face, 1) <= 3
					continue
				end
				sf = size(face, 1)
				dd = zeros(Int, sf)
				timeStamp += 1
				foreach(i -> onFace[face[i]] = timeStamp, 1 : sf)
				for i = 1 : sf
					for x in g[face[i]]
						if onFace[x] == timeStamp
							dd[i] += 1
						end
					end
				end
				mdv = 1
				for k = 2 : sf
					if dd[k] < dd[mdv]
						mdv = k
					end
				end
				for k = 1 : sf
					push!(face, face[k])
				end
				for k = mdv+2 : mdv+sf-2
					addEdge(face[mdv], face[k])
				end
			end
		end
	end

	return EdgeList

end
