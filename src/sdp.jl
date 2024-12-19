mutable struct cosmo_para
    eps_abs::Float64
    eps_rel::Float64
    max_iter::Int64
end

cosmo_para() = cosmo_para(1e-5, 1e-5, 1e4)

function beta_number(supp, coe, n, edges, d, s; constraints=nothing, QUIET=false, solver="Mosek", cosmo_setting=cosmo_para())
    println("********************************** BetaNumber **********************************")
    println("BetaNumber is launching...")
    # mbasis = get_mbasis(n, d)
    # tbasis = get_mbasis(n, min(n, 2d))
    mbasis = get_mbasis(n, n)
    tbasis = get_mbasis(n, n)
    sort!(tbasis)
    ind = [bfind(tbasis, length(tbasis), mbasis[i]) for i = 1:length(mbasis)]
    ind0 = [bfind(tbasis, length(tbasis), [i]) for i = 1:n]
    # temp = [[UInt8[1]]; get_sbasis(Vector(2:sum(length.(mbasis).<=s)), d)[2:end]]
    temp = [[UInt8[1]]; get_tbasis(Vector(2:sum(length.(mbasis).<=s)), d)]
    temp = temp[[sum(length.(mbasis[temp[i]])) for i = 1:length(temp)] .<= d]
    basis = Vector{Vector{Vector{Vector{UInt8}}}}(undef, 2)
    basis[1] = Vector{Vector{UInt8}}[]
    basis[2] = Vector{Vector{UInt8}}[]
    for i = 1:length(temp), j = 1:length(mbasis)
        l = sum(length.(mbasis[temp[i]])) + length(mbasis[j])
        w = vcat(mbasis[temp[i]]...)
        if l <= d && reduce(w) == mbasis[j] && sum(length.(mbasis[temp[i]])) != 6 && (sum(length.(mbasis[temp[i]])) == length(reduce(w)) || l < 6) && length(temp[i]) != 3 && (l != 4 || length(mbasis[j]) != 2)
        # if 2*length(mbasis[j]) <= d && reduce(w) == mbasis[j]
            # if iseven(l)
                push!(basis[1], [temp[i], [j]])
            # else
            #     push!(basis[2], [temp[i], [j]])
            # end
        end
    end
    # for j = 1:length(mbasis)
    #     if 2*length(mbasis[j]) > d
    #         push!(basis[1], [mbasis[j].+1, [j]])
    #     end
    # end
    tsupp = Vector{UInt16}[]
    lb = length.(basis)
    println(lb)
    for k = 1:1, i = 1:lb[k], j = i:lb[k]
        word,v1 = reduce(mbasis[basis[k][i][2][1]], mbasis[basis[k][j][2][1]], edges)
        v2 = reduce(mbasis[basis[k][j][2][1]], mbasis[basis[k][i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[k][i][1] == [1] && basis[k][j][1] == [1]
                word = [1]
            elseif a == 1 && basis[k][i][1] == [1]
                word = ind[basis[k][j][1]]
            elseif a == 1
                word = sort([ind[basis[k][i][1]]; ind[basis[k][j][1]]])
            elseif basis[k][i][1] == [1] && basis[k][j][1] == [1]
                word = [a]
            elseif basis[k][i][1] == [1]
                word = sort([ind[basis[k][j][1]]; a])
            else
                word = sort([ind[basis[k][i][1]]; ind[basis[k][j][1]]; a])
            end
            push!(tsupp, word)
        end
    end
    unique!(tsupp)
    sort!(tsupp)
    ltsupp = length(tsupp)
    if solver == "COSMO"
        model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => cosmo_setting.eps_abs, "eps_rel" => cosmo_setting.eps_rel, "max_iter" => cosmo_setting.max_iter))
    else
        model = Model(optimizer_with_attributes(Mosek.Optimizer))
    end
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    for k = 1:1
        pos = @variable(model, [1:lb[k], 1:lb[k]], PSD)
        for i = 1:lb[k], j = i:lb[k]
        word,v1 = reduce(mbasis[basis[k][i][2][1]], mbasis[basis[k][j][2][1]], edges)
        v2 = reduce(mbasis[basis[k][j][2][1]], mbasis[basis[k][i][2][1]], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            if a == 1 && basis[k][i][1] == [1] && basis[k][j][1] == [1]
                word = [1]
            elseif a == 1 && basis[k][i][1] == [1]
                word = ind[basis[k][j][1]]
            elseif a == 1
                word = sort([ind[basis[k][i][1]]; ind[basis[k][j][1]]])
            elseif basis[k][i][1] == [1] && basis[k][j][1] == [1]
                word = [a]
            elseif basis[k][i][1] == [1]
                word = sort([ind[basis[k][j][1]]; a])
            else
                word = sort([ind[basis[k][i][1]]; ind[basis[k][j][1]]; a])
            end
            loc = bfind(tsupp, ltsupp, word)
            if i == j
                @inbounds add_to_expression!(cons[loc], v1, pos[i,i])
            else
                @inbounds add_to_expression!(cons[loc], 2*v1, pos[i,j])
            end
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
    if constraints !== nothing
        z = @variable(model, [1:length(constraints)], lower_bound=0)
        for i = 1:length(constraints)
            cons[1] += constraints[i][2]*z[i]
            for j = 1:length(supp)
                @inbounds add_to_expression!(cons[Locb[j]], -constraints[i][1][j], z[i])
            end
        end
    end
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

function reduce(w)
    w = sort(w)
    i = 1
    while i < length(w)
        if w[i] == w[i+1]
            deleteat!(w, i)
            deleteat!(w, i)
        else
            i += 1
        end
    end
    return w
end