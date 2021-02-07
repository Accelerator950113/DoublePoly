include("Graph.jl")
include("Compress.jl")
include("PlanarPolyDraw.jl")
include("MakePicture.jl")
include("AddLine.jl")

using Printf

function DrawNetwork(ags1, ags2)
	lg = open("../logs/DrawLog.txt", "a")

	# read network data
	filePath = "../data/"
	fileName = string(filePath, ags1, ".txt")
	G0 = readGraph(fileName, "unweighted")
	println(lg, "Network ", ags1, " has ", G0.n, " nodes and ", G0.m, " edges.")

	# get compressed network
	buf = split(ags2, ',')
	selectedNodes = (size(buf, 1) > 0) ? parse(Int, buf[1]) : 20
	selectedEdges = (size(buf, 1) > 1) ? parse(Int, buf[2]) : 20
	selectMode = (size(buf, 1) > 2) ? buf[3] : "ne"
	addMode = (size(buf, 1) > 3) ? buf[4] : "KeepOne"
	G, NodeC, maxNodeC = Compress(G0, selectedNodes, selectedEdges, selectMode, addMode)
	nci, eci, ci = CentralityInformation(G, NodeC, maxNodeC)
	println(lg, "Compress Mode : ", selectedNodes, " ", selectedEdges, " ", selectMode, " ", addMode)
	println(lg, "Compressed Graph has ", G.n, " nodes and ", G.m, " edges.")
	println(lg, "Node Centrality Information : ", @sprintf("%.2f", nci*100), " %")
	println(lg, "Edge Centrality Information : ", @sprintf("%.2f", eci*100), " %")
	println(lg, "Centrality Information : ", @sprintf("%.2f", ci*100), " %")

	# draw a planar subgraph of G by using polyline drawing
	X, Y, PATH, deletedEdgeList, L, H, preservedEdgeC, deletedEdgeC = PlanarPolyDraw(G)
	println(lg, "The planar subgraph has ", size(PATH, 1), " edges.")
	println(lg, "Planar Subgraph ECI : ", @sprintf("%.2f", (preservedEdgeC/eci)*100), " %")

	# try to add some deleted edges
	X, Y, PATH, PATH2, L, H, dropedEdges = AddLine(X, Y, PATH, deletedEdgeList, L, H)
	println(lg, "Total selected edges : ", size(PATH, 1)+size(PATH2, 1)-dropedEdges, "(+", size(deletedEdgeList, 1)-dropedEdges, ")")
	println(lg, "Dropped edges : ", dropedEdges)

	nodeInfo, edgeInfo, totalInfo = CalculatePreservedCentralityInformation(G, NodeC, maxNodeC, preservedEdgeC, deletedEdgeC, PATH2)
	println(lg, "Preserved ECI : ", @sprintf("%.2f", (edgeInfo/eci)*100), " %")
	println(lg)

	# generate picture
	generatePicture(X, Y, G.V, PATH, PATH2, NodeC, ags1, max(L, H))
end
