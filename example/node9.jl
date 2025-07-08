using BetaNumber
using PyCall

graph = include("D:/project/betanumber/node9_4_hard.txt")

io1 = open("D:/project/betanumber/node9_5.txt", "w")
io2 = open("D:/project/betanumber/node9_5_hard.txt", "w")
nx = pyimport("networkx")
n = 9
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
    println(k)
    edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
    sort!(edges)
    coe = graph[k][2]
    ub = beta_number(supp, coe, n, edges, order=3, QUIET=true)[1]
    write(io1, "[$k, $ub]\n")
    lb = graph[k][3]
    if ub - lb > 1e-5
        g = graph[k][1]
        write(io2, "[\"$g\", $coe, $lb, $ub],\n")
    end
end
close(io1)
close(io2)
