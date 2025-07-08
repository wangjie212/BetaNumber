using Pkg
cd("/home/jwang/BetaNumber")
Pkg.activate(".")
using BetaNumber
using PyCall
using Graphs

nx = pyimport("networkx")
graph = include("/home/jwang/BetaNumber/example/totest9")

io1 = open("/home/jwang/BetaNumber/example/node9_huber.txt", "w")
io2 = open("/home/jwang/BetaNumber/example/node9_wang.txt", "w")
io3 = open("/home/jwang/BetaNumber/example/node9_huber_hardcase.txt", "w")
io4 = open("/home/jwang/BetaNumber/example/node9_wang_hardcase.txt", "w")
io5 = open("/home/jwang/BetaNumber/example/node9_hardcase.txt", "w")
n = 9
supp = [[i;i] for i = 1:n]
for k = 1:length(graph)
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
