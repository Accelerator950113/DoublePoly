include("Graph.jl")

function generatePicture(X, Y, labels, PATH, PATH2, NodeC, na, maxd)
	n = size(X, 1)
	m = size(PATH, 1)
	out1 = open("nodeInfo.txt", "w")
    out2 = open("edgeInfo.txt", "w")
	outlabel = open("labelInfo.txt", "w")
	outColor = open("colorInfo.txt", "w")
	outSetting = open("DrawSettings.txt", "w")
	println(outSetting, 5*maxd)
	close(outSetting)
	println(out1, n)
    for i = 1 : n
        println(out1, X[i])
        println(out1, Y[i])
    end
    close(out1)
	for i = 1 : n
		println(outColor, NodeC[labels[i]])
	end
	close(outColor)
	println(out2, m)
	for i = 1 : m
		println(out2, size(PATH[i], 1))
		for j = 1 : size(PATH[i], 1)
			println(out2, PATH[i][j][1])
			println(out2, PATH[i][j][2])
		end
	end
	close(out2)
	for i = 1 : n
		println(outlabel, labels[i])
	end
	close(outlabel)
	outSubline = open("SublineInfo.txt", "w")
	numPath2 = 0
	for i = 1 : size(PATH2, 1)
		if size(PATH2[i], 1) > 0
			numPath2 += 1
		end
	end
	println(outSubline, numPath2)
	for i = 1 : size(PATH2, 1)
		if size(PATH2[i], 1) == 0
			continue
		end
		println(outSubline, size(PATH2[i], 1))
		for j = 1 : size(PATH2[i], 1)
			println(outSubline, PATH2[i][j][1])
			println(outSubline, PATH2[i][j][2])
		end
	end
	close(outSubline)
	run(Cmd(["mpost", "-numbersystem=double", "DrawGraph.mp"]))
	run(Cmd(["mv", "DrawGraph-1.eps", string("../pictures/", na, ".eps")]))
	rm("nodeInfo.txt")
	rm("edgeInfo.txt")
	rm("labelInfo.txt")
	rm("colorInfo.txt")
	rm("DrawSettings.txt")
	rm("SublineInfo.txt")
	rm("DrawGraph.log")
end
