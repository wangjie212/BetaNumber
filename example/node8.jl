using BetaNumber
using PyCall
using Graphs

nx = pyimport("networkx")
graph = include("D:/project/betanumber/totest8")

io1 = open("D:/project/betanumber/node8_huber.txt", "w")
io2 = open("D:/project/betanumber/node8_wang.txt", "w")
io3 = open("D:/project/betanumber/node8_huber_hardcase.txt", "w")
io4 = open("D:/project/betanumber/node8_wang_hardcase.txt", "w")
io5 = open("D:/project/betanumber/node8_hardcase.txt", "w")
n = 8
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
    println(k)
    edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
    sort!(edges)
    coe = graph[k][2]
    ub1 = beta_number(supp, coe, n, edges, type="huber", QUIET=true)[1]
    write(io1, "[$k, $ub1]\n")
    ub2 = beta_number(supp, coe, n, edges, QUIET=true)[1]
    write(io2, "[$k, $ub2]\n")
    lb = max(graph[k][3], graph[k][4])
    if ub1 - lb > 1e-5
        g = graph[k][1]
        write(io3, "[\"$g\", $coe, $lb, $ub1]\n")
    end
    if ub2 - lb > 1e-5
        g = graph[k][1]
        write(io4, "[\"$g\", $coe, $lb, $ub2]\n")
    end
    ub = min(ub1, ub2)
    if ub - lb > 1e-5
        g = graph[k][1]
        write(io5, "[\"$g\", $coe, $lb, $ub]\n")
    end
end
close(io1)
close(io2)
close(io3)
close(io4)
close(io5)
