function beta_number(supp, coe, n, edges; QUIET=false)
    println("********************************** BetaNumber **********************************")
    println("BetaNumber is launching...")
    mbasis = get_mbasis(n, 3)
    tbasis = get_mbasis(n, 6)
    sort!(tbasis)
    ind = [bfind(tbasis, length(tbasis), mbasis[i]) for i = 1:length(mbasis)]
    ind0 = [bfind(tbasis, length(tbasis), [i]) for i = 1:n]
    temp = [[UInt8[1]]; get_tbasis(Vector(2:sum(length.(mbasis).<=2)), 2)]
    temp = temp[[sum(length.(mbasis[temp[i]])) for i = 1:length(temp)] .<= 3]
    basis = Vector{Vector{UInt8}}[]
    for i = 1:length(temp), j = 1:length(mbasis)
        l = sum(length.(mbasis[temp[i]])) + length(mbasis[j])
        w = vcat(mbasis[temp[i]]...)
        if l <= 6 && reduce(w) == mbasis[j] && sum(length.(mbasis[temp[i]])) != 6 && (sum(length.(mbasis[temp[i]])) == length(reduce(w)) || l < 6) && length(temp[i]) != 3 && (l != 4 || length(mbasis[j]) != 2)
            push!(basis, [temp[i], [j]])
        end
    end
    tsupp = Vector{UInt16}[]
    lb = length(basis)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(mbasis[basis[i][2][1]], mbasis[basis[j][2][1]], edges)
        v2 = reduce(mbasis[basis[j][2][1]], mbasis[basis[i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[i][1] == [1] && basis[j][1] == [1]
                word = [1]
            elseif a == 1 && basis[i][1] == [1]
                word = ind[basis[j][1]]
            elseif a == 1
                word = sort([ind[basis[i][1]]; ind[basis[j][1]]])
            elseif basis[i][1] == [1] && basis[j][1] == [1]
                word = [a]
            elseif basis[i][1] == [1]
                word = sort([ind[basis[j][1]]; a])
            else
                word = sort([ind[basis[i][1]]; ind[basis[j][1]]; a])
            end
            push!(tsupp, word)
        end
    end
    unique!(tsupp)
    sort!(tsupp)
    ltsupp = length(tsupp)
    println("The SDP size: [$lb, $ltsupp]")
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = @variable(model, [1:lb, 1:lb], PSD)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(mbasis[basis[i][2][1]], mbasis[basis[j][2][1]], edges)
        v2 = reduce(mbasis[basis[j][2][1]], mbasis[basis[i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[i][1] == [1] && basis[j][1] == [1]
                word = [1]
            elseif a == 1 && basis[i][1] == [1]
                word = ind[basis[j][1]]
            elseif a == 1
                word = sort([ind[basis[i][1]]; ind[basis[j][1]]])
            elseif basis[i][1] == [1] && basis[j][1] == [1]
                word = [a]
            elseif basis[i][1] == [1]
                word = sort([ind[basis[j][1]]; a])
            else
                word = sort([ind[basis[i][1]]; ind[basis[j][1]]; a])
            end
            loc = bfind(tsupp, ltsupp, word)
            if i == j
                @inbounds add_to_expression!(cons[loc], v1, pos[i,i])
            else
                @inbounds add_to_expression!(cons[loc], 2*v1, pos[i,j])
            end
        end
    end
    bc = zeros(ltsupp)
    for i = 1:length(supp)
        Locb = bfind(tsupp, ltsupp, ind0[supp[i]])
        bc[Locb] = -coe[i]
    end
    @variable(model, up)
    cons[1] -= up
    Locb = [bfind(tsupp, ltsupp, ind0[item]) for item in supp]
    @constraint(model, con[i=1:ltsupp], cons[i]==bc[i])
    @objective(model, Min, up)
    if QUIET == false
        println("Solving the SDP...")
    end
    optimize!(model)
    status = termination_status(model)
    objv = objective_value(model)
    if status != MOI.OPTIMAL
       println("termination status: $status")
       status = primal_status(model)
       println("solution status: $status")
    end
    println("optimum = $objv")
    dual_var = -dual.(con)
    return objv,dual_var[Locb]
end
