function beta_number(supp, coe, n, edges; type="wang", d=n, basis=[], QUIET=false)
    println("********************************** BetaNumber **********************************")
    println("BetaNumber is launching...")
    tbasis = get_mbasis(n, n)
    sort!(tbasis)
    ind = [bfind(tbasis, length(tbasis), [i]) for i = 1:n]
    if isempty(basis)
        basis = [[UInt16[1], UInt16[1]]]
        if type == "huber"
            for i = 2:length(tbasis)
                if length(tbasis[i]) <= d
                    push!(basis, [ind[tbasis[i]], [i]])
                end
            end
        else
            for i = 1:n
                push!(basis, [[ind[i]], [ind[i]]])
            end
            for i = 1:n-1, j = i+1:n
                a = bfind(tbasis, length(tbasis), [i;j])
                push!(basis, [[ind[i]; a], [ind[j]]], [[ind[j]; a], [ind[i]]])
            end
            for i = 1:n-1, j = i+1:n, k in setdiff(Vector(1:n), [i;j])
                a = bfind(tbasis, length(tbasis), [i;j])
                b = bfind(tbasis, length(tbasis), sort([i;j;k]))
                push!(basis, [sort([a; ind[k]]), [b]])
                # c = bfind(tbasis, length(tbasis), sort([i;k]))
                # push!(basis, [sort([a; ind[j]; ind[k];]), [c]])
            end
        end
    end
    tsupp = Vector{UInt16}[]
    lb = length(basis)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(tbasis[basis[i][2][1]], tbasis[basis[j][2][1]], edges)
        v2 = reduce(tbasis[basis[j][2][1]], tbasis[basis[i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[i][1] == [1] && basis[j][1] == [1]
                word = [1]
            elseif a == 1 && basis[i][1] == [1]
                word = basis[j][1]
            elseif a == 1
                word = sort([basis[i][1]; basis[j][1]])
            elseif basis[i][1] == [1] && basis[j][1] == [1]
                word = [a]
            elseif basis[i][1] == [1]
                word = sort([basis[j][1]; a])
            else
                word = sort([basis[i][1]; basis[j][1]; a])
            end
            push!(tsupp, word)
        end
    end
    unique!(tsupp)
    sort!(tsupp)
    ltsupp = length(tsupp)
    println("The SDP size: [$lb, $ltsupp]")
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_IPAR_NUM_THREADS" =>8))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    pos = @variable(model, [1:lb, 1:lb], PSD)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(tbasis[basis[i][2][1]], tbasis[basis[j][2][1]], edges)
        v2 = reduce(tbasis[basis[j][2][1]], tbasis[basis[i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[i][1] == [1] && basis[j][1] == [1]
                word = [1]
            elseif a == 1 && basis[i][1] == [1]
                word = basis[j][1]
            elseif a == 1
                word = sort([basis[i][1]; basis[j][1]])
            elseif basis[i][1] == [1] && basis[j][1] == [1]
                word = [a]
            elseif basis[i][1] == [1]
                word = sort([basis[j][1]; a])
            else
                word = sort([basis[i][1]; basis[j][1]; a])
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
        Locb = bfind(tsupp, ltsupp, ind[supp[i]])
        bc[Locb] = -coe[i]
    end
    @variable(model, up)
    cons[1] -= up
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
    moment = zeros(lb, lb)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(tbasis[basis[i][2][1]], tbasis[basis[j][2][1]], edges)
        v2 = reduce(tbasis[basis[j][2][1]], tbasis[basis[i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[i][1] == [1] && basis[j][1] == [1]
                word = [1]
            elseif a == 1 && basis[i][1] == [1]
                word = basis[j][1]
            elseif a == 1
                word = sort([basis[i][1]; basis[j][1]])
            elseif basis[i][1] == [1] && basis[j][1] == [1]
                word = [a]
            elseif basis[i][1] == [1]
                word = sort([basis[j][1]; a])
            else
                word = sort([basis[i][1]; basis[j][1]; a])
            end
            loc = bfind(tsupp, ltsupp, word)
            moment[i,j] = moment[j,i] = dual_var[loc]
        else
            moment[i,j] = moment[j,i] = 0
        end
    end
    return objv,moment,basis
end
