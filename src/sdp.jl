mutable struct cosmo_para
    eps_abs::Float64
    eps_rel::Float64
    max_iter::Int64
end

cosmo_para() = cosmo_para(1e-5, 1e-5, 1e4)

function beta_number(supp, coe, n, edges, d; QUIET=false, solver="Mosek", cosmo_setting=cosmo_para())
    println("********************************** BetaNumber **********************************")
    println("BetaNumber is launching...")
    basis = get_mbasis(n, min(n, d))
    tbasis = get_mbasis(n, min(n, 2d))
    sort!(tbasis)
    ind = [bfind(tbasis, length(tbasis), [i]) for i = 1:n]
    lb = length(basis)
    tsupp = Vector{UInt16}[]
    for i = 1:lb, j = i:lb
        word,v1 = reduce(basis[i], basis[j], edges)
        v2 = reduce(basis[j], basis[i], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            word = a == 1 ? sort([ind[basis[i]]; ind[basis[j]]]) : sort([ind[basis[i]]; ind[basis[j]]; a])
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
    pos = @variable(model, [1:lb, 1:lb], PSD)
    for i = 1:lb, j = i:lb
        word,v1 = reduce(basis[i], basis[j], edges)
        v2 = reduce(basis[j], basis[i], edges)[2]
        if v1 == v2
            a = bfind(tbasis, length(tbasis), word)
            word = a == 1 ? sort([ind[basis[i]]; ind[basis[j]]]) : sort([ind[basis[i]]; ind[basis[j]]; a])
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
    @constraint(model, cons.==bc)
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
    return objv
end
