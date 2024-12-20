function get_mbasis(n, d)
    nbasis = [UInt8[i] for i = 1:n]
    basis = [[UInt8[]]; nbasis]
    for i = 2:d
        nbasis = lift_basis(nbasis, n)
        append!(basis, nbasis)
    end
    return basis
end

function lift_basis(basis, n)
    nbasis = Vector{UInt8}[]
    for item in basis
        if item[end] < n
            append!(nbasis, [[item; j] for j = item[end]+1:n])
        end
    end
    return nbasis
end

function bfind(A, l, a)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] ==  a
           return mid
        elseif A[mid] < a
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end

function reduce(word1, word2, edges)
    word = [reverse(word1); word2]
    val = 1
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, [i, i+1])
            i = i == 1 ? 1 : i-1
        elseif word[i] > word[i+1]
            word[i],word[i+1] = word[i+1],word[i]
            if bfind(edges, length(edges), word[i:i+1]) !== nothing
                val *= -1
            end
            i = i == 1 ? 2 : i-1
        else
            i += 1
        end
    end
    return word,val
end
