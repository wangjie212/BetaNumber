function get_mbasis(n, d)
    nbasis = [UInt8[i] for i = 1:n]
    basis = [[UInt8[]]; nbasis]
    for i = 2:d
        nbasis = lift_basis(nbasis, n)
        append!(basis, nbasis)
    end
    return basis
end

function get_tbasis(var, d)
    basis = [UInt8[i] for i in var]
    nbasis = [UInt8[i] for i in var]
    for i = 2:d
        nbasis = lift_tbasis(nbasis, var)
        append!(basis, nbasis)
    end
    return basis
end

function lift_tbasis(basis, var)
    nbasis = Vector{UInt8}[]
    for item in basis
        if item[end] < var[end]
            k = bfind(var, length(var), item[end])
            append!(nbasis, [[item; j] for j in var[k+1:end]])
        end
    end
    return nbasis
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

function get_sbasis(var, d)
    n = length(var)
    lb = binomial(n+d, d)
    basis = Vector{Vector{UInt16}}(undef, lb)
    basis[1] = UInt16[]
    i = 0
    t = 1
    while i < d+1
        t += 1
        if sum(basis[t-1]) == var[n]*i
           if i < d
               basis[t] = var[1]*ones(UInt16, i+1)
           end
           i += 1
        else
            j = bfind(var, n, basis[t-1][1])
            basis[t] = copy(basis[t-1])
            ind = findfirst(x->basis[t][x]!=var[j], 1:length(basis[t]))
            if ind === nothing
                ind = length(basis[t])+1
            end
            if j != 1
                basis[t][1:ind-2] = var[1]*ones(UInt16, ind-2)
            end
            basis[t][ind-1] = var[j+1]
        end
    end
    return basis
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
