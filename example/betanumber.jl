using BetaNumber
# using PyCall

## graph6 code to igraph
# nx = pyimport("networkx")
# ig6 = "GUZv\\{"
# edges = [[i+1,j+1] for (i,j) in nx.from_graph6_bytes(codeunits(ig6)).edges()]
# edges = [[1,2], [1,4], [1,5], [1,7], [2,3], [2,4], [2,5], [3,4], [3,6], [4,5], [4,7], [5,6], [5,7], [6,7]]
edges = [[1,2], [2,3], [3,4], [4,5], [1,5], [1,6], [2,6], [3,6], [4,6], [5,6]]
# edges = [[1, 3], [1, 5], [1, 6], [1, 7], [2, 4], [2, 5], [2, 7], [3, 5], [3, 6], [3, 7], [4, 6], [4, 7], [5, 7], [6, 7]]
sort!(edges)
n = 6
supp = [[i;i] for i = 1:n]
coe = [ones(5); 2]
# coe = ones(n)
# coe = [1, 1, 1, 1, 1, 1, 1, 1]
beta_number(supp, coe, n, edges, QUIET=false)
