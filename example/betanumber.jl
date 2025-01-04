using BetaNumber
using PyCall
using Graphs

graph = include("E:/project/betanumber/node9_hardcase.txt")

io = open("E:/project/betanumber/node9_32.txt", "w")

## graph6 code to igraph
nx = pyimport("networkx")
n = 9
supp = [[i;i] for i = 1:n]
# for k = 1:length(graph), l = 1:length(graph[k])-1
for k = 25:40
# v = val[s]
# k = Int(v[1])
# l = Int(v[2])
edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(graph[k][1])).edges()]
sort!(edges)
coe = graph[k][2]
println(k)
v,moment,basis = beta_number(supp, coe, n, edges, type="huber", QUIET=true)
G = SimpleGraph(size(moment,2))
for i = 1:size(moment,2), j = i+1:size(moment,2)
    if abs(moment[i,j]) > 1e-6
        add_edge!(G, i, j)
    end
end
blocks = connected_components(G)
b1 = deepcopy(basis[blocks[1]])
v,moment,basis = beta_number(supp, coe, n, edges, QUIET=true)
G = SimpleGraph(size(moment,2))
for i = 1:size(moment,2), j = i+1:size(moment,2)
    if abs(moment[i,j]) > 1e-6
        add_edge!(G, i, j)
    end
end
blocks = connected_components(G)
b2 = deepcopy(basis[blocks[1]])
ba = unique([b1; b2])
v,moment,basis = beta_number(supp, coe, n, edges, basis=ba, QUIET=true)
write(io, "[$k, $v]\n")
end
close(io)

# io = open("E:/project/betanumber/moment.txt", "w")
# write(io, "$moment")
# close(io)

ind = Int[]
for (i,v) in enumerate(val)
    if abs(max(graph[Int(v[1])][Int(v[2])+1][2], graph[Int(v[1])][Int(v[2])+1][3]) - v[3]) > 1e-5
        push!(ind, i)
    end 
end

val = val[ind]

io = open("E:/project/betanumber/node9_hardcase.txt", "w")
for s = 1:69
v = val[s]
k = Int(v[1])
l = Int(v[2])
# println(v[3])
g = graph[k][1]
coe = graph[k][l+1][1]
lb = max(graph[Int(v[1])][Int(v[2])+1][2], graph[Int(v[1])][Int(v[2])+1][3])
ub = min(v[3], val1[s][3])
# println(graph[k][l+1][2])
# println(graph[k][l+1][3])
if ub - lb > 1e-5
write(io, "[\"$g\", $coe, $lb, $ub]\n")
end
end
close(io)