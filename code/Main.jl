include("Draw.jl")

tl = open("../logs/timeLog.txt", "a")
T = time()
ags2 = (size(ARGS, 1) < 2) ? "20,20" : ARGS[2]
DrawNetwork(ARGS[1], ags2)
T = time() - T
println(tl, ARGS[1], " ", T)
