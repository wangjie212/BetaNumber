using BetaNumber
# using PyCall

## graph6 code to igraph
# nx = pyimport("networkx")
# ig6 = "GUzrvW"
# edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(ig6)).edges()]
edges = [[1,2], [1,4], [1,5], [1,7], [2,3], [2,4], [2,5], [3,4], [3,6], [4,5], [4,7], [5,6], [5,7], [6,7]]
sort!(edges)
n = 7
supp = [[i;i] for i = 1:n]
coe = ones(n)
d = 3
@time opt = beta_number(supp, coe, n, edges, d)
