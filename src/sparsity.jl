function get_graph(basis, tsupp, tbasis, edges)
    lb = length(basis)
    G = SimpleGraph(lb)
    ltsupp = length(tsupp)
    for i = 1:lb, j = i+1:lb
        word = sadd(basis[i], basis[j], tbasis, edges)[1]
        if word !== nothing && bfind(tsupp, ltsupp, word) !== nothing
            add_edge!(G, i, j)
        end
    end
    return G
end

function get_blocks(basis, tsupp, tbasis, edges; QUIET=true)
    G = get_graph(basis, tsupp, tbasis, edges)
    blocks = connected_components(G)
    blocksize = length.(blocks)
    cl = length(blocksize)
    if QUIET == false
        sb = sort(unique(blocksize), rev=true)
        numb = [sum(blocksize.== i) for i in sb]
        println("-----------------------------------------------------------------------------")
        println("The sizes of PSD blocks:\n$sb\n$numb")
        println("-----------------------------------------------------------------------------")
    end
    return blocks,cl,blocksize
end
