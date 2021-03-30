# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]

function star()
    #CREATE MODEL
    feas = Model(get_optimizer())
    set_optimizer_attribute(feas,"MIPGap",0.001)

    #CREATE DELIVERY VARIABLES AND CONSTRAINTS
    R = Dict()

    for k in K(), t in T()
        u = @variable(feas, [i = K(k).cover, j = K(k).cover])
        v = @variable(feas, [j = K(k).cover])
        y = @variable(feas, [i = K(k).cover, j = K(k).cover], Int)
        z = @variable(feas, [j = K(k).cover], Int)

        @constraint(feas, [j = K(k).cover], z[j] >= 0)
        @constraint(feas, [j = K(k).cover], v[j] >= 0)

        @constraint(feas, [i = K(k).cover, j = K(k).cover], u[i,j] >= 0)
        @constraint(feas, [i = K(k).cover, j = K(k).cover], y[i,j] >= 0)

        @constraint(feas, [j = K(k).cover], z[j] <= K(k).BP[j])

        @constraint(feas, [j = K(k).cover], sum(u[i,j] for i in K(k).cover) - v[j] == 0)
        @constraint(feas, [i = K(k).cover, j = K(k).cover], u[i,j] <= K(k).Q * y[i,j])
        @constraint(feas, [j = K(k).cover], v[j] <= K(k).Q * z[j])

        R[(k,t)] = aggr(u,v,y,z)
    end

    #CREATE INVENTORY LEVEL VARIABLE AND CONSTRAINTS
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())])

    @constraint(feas, λ[i = V(), t = T()],
        I[i,t-1] + sum(sum(R[(k,t)].u[i,j] for j in K(k).cover) -
        R[(k,t)].v[i] for k in passes(i)) == d(i,t) + I[i,t]
    )

    @constraint(feas, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    )

    @constraint(feas, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    )

    #ADD APPROX OBJECTIVE VALUE
    @objective(feas, Min,
        sum(V(i).h * I[i,t] for i in V(), t in T()) +
        sum(2 * K(k).vx * dist(i,j) * R[(k,t)].y[i,j]
            for k in K(), i in K(k).cover, j in K(k).cover, t in T()
        ) +
        sum((K(k).vl * dist(i,j) + K(k).fd) * R[(k,t)].u[i,j]
            for k in K(), i in K(k).cover, j in K(k).cover, t in T()
        ) +
        sum(K(k).fp * R[(k,t)].z[j] for k in K(), j in K(k).cover, t in T())
    )

    optimize!(feas)

    return feas,R
end

function dis(R::Dict,j_dis::Int64,k_dis::Int64,t_dis::Int64)
    #CREATE DISAGGREGATED INDEX
    dis_idx = [j_dis]
    for i in K(k_dis).cover
        if value.(R[(k_dis,t_dis)].y[i,j_dis]) > 0
            push!(dis_idx,i)
        end
    end

    #CREATE AGGREGATE GUIDE (FOR CONSTRAINTS)
    v_guide = Dict{Int64,Int64}(dis_idx .=> [0 for i in dis_idx])
    v_guide[j_dis] = round(value(R[(k_dis,t_dis)].v[j_dis]))
    u_guide = Dict{Int64,Int64}(
        dis_idx .=> [round(value(R[(k_dis,t_dis)].u[i,j_dis])) for i in dis_idx]
    )
    y_guide = Dict{Int64,Int64}(
        dis_idx .=> [round(value(R[(k_dis,t_dis)].y[i,j_dis])) for i in dis_idx]
    )

    #CREATE MANIFESTS OF DISAGGREGATION
    M = collect(1:Int64(value(R[(k_dis,t_dis)].z[j_dis])))

    #CREATE MODEL
    feas = Model(get_optimizer())

    #CREATE DELIVERY VARIABLES AND CONSTRAINTS
    R = Dict()

    for m in M
        z = @variable(feas, [j = dis_idx], Bin)
        y = @variable(feas, [i = dis_idx], Bin)
        x = @variable(feas, [i = dis_idx, j = dis_idx], Bin)
        v = @variable(feas, [j = dis_idx], Int)
        u = @variable(feas, [i = dis_idx], Int)
        w = @variable(feas, [i = dis_idx, j = dis_idx], Int)

        @constraint(feas, [j = dis_idx], v[j] >= 0)
        @constraint(feas, [i = dis_idx], u[i] >= 0)
        @constraint(feas, [i = dis_idx, j = dis_idx], w[i,j] >= 0)

        @constraint(feas, sum(u[i] for i in dis_idx) - v[j_dis] == 0)
        @constraint(feas, [i = dis_idx],
            sum(w[j,i] for j in dis_idx) - sum(w[i,j] for j in dis_idx) == u[i] - v[i]
        )
        @constraint(feas, [i = dis_idx],
            sum(x[j,i] for j in dis_idx) - sum(x[i,j] for j in dis_idx) == 0
        )

        @constraint(feas, [i = dis_idx], u[i] <= K(k_dis).Q * y[i])
        @constraint(feas, v[j_dis] <= K(k_dis).Q * z[j_dis])
        @constraint(feas, [i = dis_idx, j = dis_idx], w[i,j] <= K(k_dis).Q * x[i,j])

        R[m] = disaggr(u,v,w,y,z,x)
    end

    #AGGREGATE LIMIT CONSTRAINTS
    @constraint(feas, [j = dis_idx], sum(R[m].v[j] for m in M) == v_guide[j])
    @constraint(feas, [i = dis_idx], sum(R[m].u[i] for m in M) == u_guide[i])
    @constraint(feas, [i = dis_idx], sum(R[m].y[i] for m in M) == y_guide[i])

    #OBJECTIVE FUNCTION
    @objective(feas, Min,
        sum(K(k_dis).vx * dist(i,j) * R[m].x[i,j]
            for i in dis_idx, j in dis_idx, m in M
        ) +
        sum(K(k_dis).vl * dist(i,j) * R[m].w[i,j]
            for i in dis_idx, j in dis_idx, m in M
        ) +
        sum(K(k_dis).fd * R[m].u[i]
            for i in dis_idx, m in M
        ) +
        sum(K(k_dis).fp * R[m].z[j]
            for j in dis_idx, m in M
        )
    )

    optimize!(feas)

    return feas,R
end
