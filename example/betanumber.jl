using BetaNumber
using PyCall

nx = pyimport("networkx")
graph = include("D:/project/betanumber/node8_hardcase.txt")
io = open("D:/project/betanumber/node8_hardcase2.txt", "w")

n = 8
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
println(k)
edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
sort!(edges)
coe = graph[k][2]
v = beta_number(supp, coe, n, edges, order=3, QUIET=false)[1]
write(io, "[$k, $v]\n")
end
close(io)

graph1 = include("D:/project/betanumber/node9_hardcase.txt")
graph2 = include("D:/project/betanumber/node9_hardcase2.txt")
io = open("D:/project/betanumber/node9_hardcase3.txt", "w")
for k = 1:length(graph1)
    if graph2[k][2] - graph1[k][3] > 1e-5
        g = graph1[k][1]
        coe = graph1[k][2]
        lb = graph1[k][3]
        ub = graph2[k][2]
        write(io, "[\"$g\", $coe, $lb, $ub],\n")
    end
end
close(io)

graph = include("D:/project/betanumber/node9_7_hard.txt")
io = open("D:/project/betanumber/node9_7_hard1.txt", "w")
for k = 1:length(graph)
    if graph[k][4] - graph[k][3] > 1e-3
        g = graph[k][1]
        coe = graph[k][2]
        lb = graph[k][3]
        ub = graph[k][4]
        write(io, "[\"$g\", $coe, $lb, $ub],\n")
    end
end
close(io)
