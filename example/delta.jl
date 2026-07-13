using BetaNumber
using PyCall

nx = pyimport("networkx")
graph = include("D:/project/betanumber/delta/test8.txt")

# io1 = open("D:/project/betanumber/delta/ndelta8_2.txt", "w")
# io2 = open("D:/project/betanumber/delta/ndelta8_2_hard.txt", "w")

n = 8
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
    println(k)
    edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
    sort!(edges)
    a = delta_number(supp, n, edges, order=3, QUIET=false)[1]
    ub = 1/a
    val = 1/graph[k][2]
    write(io1, "[$k, $val, $ub]\n")
    if abs(ub - val) > 1e-5
        g = graph[k][1]
        write(io2, "[$k, $val, $ub],\n")
    end
end
close(io1)
close(io2)
