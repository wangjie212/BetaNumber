function beta_number(supp, coe, n, edges; order=2, QUIET=false)
    println("********************************** BetaNumber **********************************")
    println("BetaNumber is launching...")
    tbasis = get_mbasis(n, n)
    sort!(tbasis)
    ind = [bfind(tbasis, length(tbasis), [i]) for i = 1:n]
    basis = [[UInt16[1], UInt16[1]]]
    for i = 1:n
        push!(basis, [[ind[i]], [ind[i]]])
    end
    if order > 1
        for i = 1:n-1, j = i+1:n
            a = bfind(tbasis, length(tbasis), [i;j])
            push!(basis, [sort([ind[i]; a]), [ind[j]]], [sort([ind[j]; a]), [ind[i]]])
            push!(basis, [[ind[i]; ind[j]], [a]])
        end
    end
    if order > 2
        for i = 1:n-1, j = i+1:n, k in setdiff(Vector(1:n), [i;j])
            a = bfind(tbasis, length(tbasis), [i;j])
            b = bfind(tbasis, length(tbasis), sort([i;j;k]))
            push!(basis, [sort([a; ind[k]]), [b]])
            c = bfind(tbasis, length(tbasis), sort([i;k]))
            push!(basis, [sort([a; ind[j]; ind[k]]), [c]])
            push!(basis, [sort([a; ind[i]; ind[j]; ind[k]]), [ind[k]]])
        end
        for i = 1:n-2, j = i+1:n-1, k = j+1:n
            a = bfind(tbasis, length(tbasis), [i;j;k])
            push!(basis, [[ind[i]; ind[j]; ind[k]], [a]])
        end
    end
    if order > 3
        for i = 1:n-1, j = i+1:n, k in setdiff(Vector(1:n-1), [i;j]), l in setdiff(Vector(k+1:n), [i;j])
            a = bfind(tbasis, length(tbasis), [i;j])
            b = bfind(tbasis, length(tbasis), sort([i;j;k;l]))
            push!(basis, [sort([a; ind[k]; ind[l]]), [b]])
            # c = bfind(tbasis, length(tbasis), sort([j;k;l]))
            # push!(basis, [sort([a; ind[i]; ind[k]; ind[l]]), [c]])
            # c = bfind(tbasis, length(tbasis), sort([i;k;l]))
            # push!(basis, [sort([a; ind[j]; ind[k]; ind[l]]), [c]])
        end
        # for i = 1:n-3, j = i+1:n-2, k = j+1:n-1, l = k+1:n
        #     a = bfind(tbasis, length(tbasis), [i;j;k;l])
        #     push!(basis, [[ind[i]; ind[j]; ind[k]; ind[l]], [a]])
        # end
    end
    tsupp = [sadd(item, item, tbasis, edges)[1] for item in basis]
    unique!(tsupp)
    sort!(tsupp)
    blocks,cl,blocksize = get_blocks(basis, tsupp, tbasis, edges, QUIET=QUIET)
    tsupp = Vector{UInt16}[]
    for k = 1:cl, i = 1:blocksize[k], j = i:blocksize[k]
        word = sadd(basis[blocks[k][i]], basis[blocks[k][j]], tbasis, edges)[1]
        if word !== nothing
            push!(tsupp, word)
        end
    end
    unique!(tsupp)
    sort!(tsupp)
    ltsupp = length(tsupp)
    println("There are $ltsupp affine constraints.")
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_IPAR_NUM_THREADS" =>4))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    for k = 1:cl
        if blocksize[k] == 1
            pos = @variable(model, lower_bound=0)
            word = sadd(basis[blocks[k][1]], basis[blocks[k][1]], tbasis, edges)[1]
            loc = bfind(tsupp, ltsupp, word)
            @inbounds add_to_expression!(cons[loc], pos)
        else
            pos = @variable(model, [1:blocksize[k], 1:blocksize[k]], PSD)
            for i = 1:blocksize[k], j = i:blocksize[k]
                word,v = sadd(basis[blocks[k][i]], basis[blocks[k][j]], tbasis, edges)
                if word !== nothing
                    loc = bfind(tsupp, ltsupp, word)
                    if i == j
                        @inbounds add_to_expression!(cons[loc], pos[i,i])
                    else
                        @inbounds add_to_expression!(cons[loc], 2*v, pos[i,j])
                    end
                end
            end
        end
    end
    for i = 1:length(supp)
        Locb = bfind(tsupp, ltsupp, ind[supp[i]])
        cons[Locb] += coe[i]
    end
    @variable(model, up)
    cons[1] -= up
    @constraint(model, con, cons==zeros(ltsupp))
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
    dual_var = -dual(con)
    moment = Vector{Matrix{Float64}}(undef, cl)
    for k = 1:cl
        moment[k] = zeros(blocksize[k], blocksize[k])
        for i = 1:blocksize[k], j = i:blocksize[k]
            word,v = sadd(basis[blocks[k][i]], basis[blocks[k][j]], tbasis, edges)
            if word !== nothing
                loc = bfind(tsupp, ltsupp, word)
                moment[k][i,j] = moment[k][j,i] = v*dual_var[loc]
            else
                moment[k][i,j] = moment[k][j,i] = 0
            end
        end
    end
    return objv,moment,basis
end
