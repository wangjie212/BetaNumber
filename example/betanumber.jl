using BetaNumber
using PyCall

nx = pyimport("networkx")
graph = include("D:/project/betanumber/node9_hard12.txt")
io = open("D:/project/betanumber/node8_hardcase2.txt", "w")

n = 9
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
println(k)
edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
sort!(edges)
println(edges)
coe = graph[k][2]
v = beta_number(supp, coe, n, edges, order=2, QUIET=false)[1]
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

graph = include("D:/project/betanumber/node9_hard5.txt")
io = open("D:/project/betanumber/node9_hard12.txt", "w")
for k = 1:length(graph)
    if graph[k][4] - graph[k][3] > 4.7e-4
        g = graph[k][1]
        coe = graph[k][2]
        lb = graph[k][3]
        ub = graph[k][4]
        write(io, "[\"$g\", $coe, $lb, $ub],\n")
    end
end
close(io)


graph = include("D:/project/betanumber/delta/test8.txt")

io1 = open("/home/jwang/Programs/delta/ndelta8_1.txt", "w")
io2 = open("/home/jwang/Programs/delta/ndelta8_1_hard.txt", "w")

n = 8
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
    println(k)
    edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
    sort!(edges)
    ub = delta_number(supp, n, edges, order=3, QUIET=false)[1]
    val = 1/graph[k][2]
    write(io1, "[$k, $val, $ub]\n")
    if abs(ub - val) > 1e-5
        g = graph[k][1]
        write(io2, "[$k, $val, $ub],\n")
    end
    println(1/v)
end
close(io1)
close(io2)
