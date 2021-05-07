# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================

passes(i) = [k for k in K() if (i in K(k).cover)]

function star()
    #CREATE MODEL
    feas = Model(get_optimizer())
    #set_optimizer_attribute(feas,"MIPGap",0.005)
    #set_optimizer_attribute(feas,"Presolve",2)
    #set_optimizer_attribute(feas,"MIPFocus",3)

    #CREATE DELIVERY VARIABLES AND CONSTRAINTS
    R = Dict()

    for k in K(), t in T()
        u = @variable(feas, [i = K(k).cover])#, Int)
        v = @variable(feas, [j = K(k).cover])#, Int)
        y = @variable(feas, [i = K(k).cover])#, Int)
        z = @variable(feas, [j = K(k).cover])#, Int)
        o = @variable(feas, [i = K(k).cover, j = K(k).cover])#, Int)
        p = @variable(feas, [i = K(k).cover, j = K(k).cover])#, Int)

        @constraint(feas, [j = K(k).cover], z[j] >= 0)
        @constraint(feas, [j = K(k).cover], v[j] >= 0)
        @constraint(feas, [i = K(k).cover], u[i] >= 0)
        @constraint(feas, [i = K(k).cover], y[i] >= 0)
        @constraint(feas, [i = K(k).cover, j = K(k).cover], o[i,j] >= 0)
        @constraint(feas, [i = K(k).cover, j = K(k).cover], p[i,j] >= 0)

        @constraint(feas, [j = K(k).cover], z[j] <= K(k).BP[j])

        @constraint(feas,
            sum(u[i] for i in K(k).cover) - sum(v[j] for j in K(k).cover) == 0
        ) #all pickup delivered

        @constraint(feas, [j = K(k).cover],
            sum(o[i,j] for i in K(k).cover) == v[j]
        ) #all flow out of depot is pickup

        @constraint(feas, [i = K(k).cover],
            sum(o[i,j] for j in K(k).cover) == u[i]
        ) #all flow to depot is delivery

        @constraint(feas, [i = K(k).cover], u[i] <= K(k).Q * y[i])
        @constraint(feas, [j = K(k).cover], v[j] <= K(k).Q * z[j])
        @constraint(feas, [i = K(k).cover, j = K(k).cover], o[i,j] <= K(k).Q * p[i,j])

        R[(k,t)] = aggr(u,v,y,z,o,p)
    end

    #CREATE INVENTORY LEVEL VARIABLE AND CONSTRAINTS
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())])

    @constraint(feas, λ[i = V(), t = T()],
        I[i,t-1] + sum(R[(k,t)].u[i] for k in passes(i)) ==
        d(i,t) + sum(R[(k,t)].v[i] for k in passes(i)) + I[i,t]
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
        sum(2 * K(k).vx * dist(i,j) * R[(k,t)].p[i,j]
            for k in K(), i in K(k).cover, j in K(k).cover, t in T()
        ) +
        sum(K(k).vl * dist(i,j) * R[(k,t)].o[i,j]
            for k in K(), i in K(k).cover, j in K(k).cover, t in T()
        ) +
        sum(K(k).fd * R[(k,t)].u[i] for k in K(), i in K(k).cover, t in T()) +
        sum(K(k).fp * R[(k,t)].z[j] for k in K(), j in K(k).cover, t in T())
    )

    optimize!(feas)

    return feas,R
end

function generate(R::Dict,j::Int64,k_dis::Int64,t_dis::Int64)
    println("=====START DISAGGR ON j:$j,k:$k_dis,t:$t_dis=====")
    routes = Dict{Int64,disaggr}() #SOLUTION VECTOR

    source = j
    sink = []
    for i in K(k_dis).cover
        if value.(R[(k_dis,t_dis)].o[i,source]) > 0
            push!(sink,i)
        end
    end #SINK DIDAPAT DARI AGGREGATE SOLUTION + K_DIS, T_DIS, J (SOURCE)
    nodes = union(source,sink)

    demands = Dict{Int64,Int64}(
        sink .=> [round.(value.(R[(k_dis,t_dis)].o[i,source])) for i in sink]
    ) #DEMANDS NGEEXTRACT VALUE O NYA

    NV = K(k_dis).BP[source] #NV DIDAPAT DR K_DIS, SOURCE

    if isempty(sink)
        println("no delivery")
        return routes,demands,NV
    end #IF THERE ARE NO DELIVERY TO BE MADE

    for i in sink
        n = floor(demands[i]/K(k_dis).Q) #number of direct deliveries to i
        if n > 0
            println("$n direct deliveries from $source to $i with size $(K(k_dis).Q)")
            for m in 1:n
                u = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes)
                v = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes)
                w = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes,nodes)

                y = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes)
                z = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes)
                x = JuMP.Containers.DenseAxisArray{Float64}(undef,nodes,nodes)

                u .= 0
                v .= 0
                w .= 0
                y .= 0
                z .= 0
                x .= 0

                u[i] = K(k_dis).Q
                v[source] = K(k_dis).Q
                w[source,i] = K(k_dis).Q
                w[i,source] = 0

                y[i] = 1
                z[source] = 1
                x[source,i] = 1
                x[i,source] = 1

                routes[length(keys(routes))+1] = disaggr(u,v,w,y,z,x)
            end

            NV -= n #number of vehicle reduced
            demands[i] -= K(k_dis).Q * n #demand covered

            if demands[i] == 0
                println("all demand to $i from $source covered.")
            end
        end
    end #CREATE DIRECT DELIVERIES

    println("sisa demand: $demands")
    println("sisa kendaraan: $NV")

    if sum(demands[i] for i in keys(demands)) != 0
        additional_route = residual(source,k_dis,demands,NV,K(k_dis).Q) #solve res demands
        additional_route = filter(p -> value.(last(p).ind) == 1, additional_route)

        println("added $(length(additional_route)) routes")
        for r in additional_route
            for i in last(r).u.axes[1][2:end] #source not included
                demands[i] -= round.(value.(last(r).u[i]))
            end

            NV -= value.(last(r).ind)

            routes[length(keys(routes))+1] = disaggr(
                round.(value.(last(r).u)), round.(value.(last(r).v)), round.(value.(last(r).w)),
                round.(value.(last(r).y)), round.(value.(last(r).z)), round.(value.(last(r).x))
            )
        end
    end

    return routes,demands,NV
end #GENERATE (INPUT: AGGR SOL, J, K_DIS, T_DIS)

function residual(source,k,demands,NV,Q)
    res = Model(get_optimizer())
    set_silent(res)
    sink = collect(keys(filter(p -> last(p) > 0,demands)))
    nodes = union(source,sink)

    routes = Dict{Int64,NamedTuple}()

    for m in 1:NV
        u = @variable(res, [i = nodes])
        v = @variable(res, [i = nodes])
        w = @variable(res, [i = nodes, j = nodes])
        y = @variable(res, [i = nodes], binary = true)
        z = @variable(res, [i = nodes], binary = true)
        x = @variable(res, [i = nodes, j = nodes], binary = true)

        p = @variable(res, [i = nodes]) #fraction of demand
        ind = @variable(res, binary = true)

        @constraint(res, [i = nodes], u[i] >= 0)
        @constraint(res, [i = nodes], v[i] >= 0)
        @constraint(res, [i = nodes, j = nodes], w[i,j] >= 0)

        @constraint(res, sum(z[i] for i in sink) == 0) #can only start at source

        #CONSTRAINTS FROM DROR TRUEDEAU (1990) SDVRP FOR EACH VEHICLE
        @constraint(res, [i = nodes],
            sum(x[j,i] for j in nodes) - sum(x[i,j] for j in nodes) == 0
        ) #traversal

        @constraint(res,
            sum(demands[i] * p[i] for i in sink) <= Q
        ) #total of demand in vehicle

        @constraint(res, [i = sink],
            p[i] >= 0
        ) #fraction >= 0

        #CONSTRAINTS FROM ORIGINAL FORMULATION OF INVENTORY ROUTING + CONNECTING FRACTION
        @constraint(res, [i = sink],
            u[i] == p[i] * demands[i]
        ) #u is fraction of demand fulfilled

        @constraint(res,
            sum(u[i] for i in nodes) - sum(v[i] for i in nodes) == 0
        ) #all pickup delivered

        @constraint(res, [i = nodes],
            sum(w[j,i] for j in nodes) - sum(w[i,j] for j in nodes) == u[i] - v[i]
        ) #vehicle load

        @constraint(res, [i = nodes], u[i] <= Q * y[i]) #delivery correlation
        @constraint(res, [i = nodes], v[i] <= Q * z[i]) #pickup correlation
        @constraint(res, [i = nodes, j = nodes], w[i,j] <= Q * x[i,j]) #travel correl

        @constraint(res, [i = nodes], y[i] + z[i] <= 1)

        @constraint(res, [i = nodes],
            sum(x[j,i] for j in nodes) == y[i] + z[i]
        ) #route in

        @constraint(res, [i = nodes],
            sum(x[i,j] for j in nodes) == y[i] + z[i]
        ) #route out

        @constraint(res, [i = nodes],
            y[i] <= ind
        )

        @constraint(res, [i = nodes],
            z[i] <= ind
        )

        routes[m] = (u = u, v = v, w = w, y = y, z = z, x = x, p = p, ind = ind)
    end

    #CONSTRAINT2 GABUNGAN FROM DROR TRUEDEAU (1990) SDVRP
    @constraint(res, [i = sink], sum(routes[m].p[i] for m in 1:NV) == 1)

    #OBJECTIVE FUNCTION ASLI INVENTORY ROUTING
    @objective(res, Min,
        sum(K(k).vx * dist(i,j) * routes[m].x[i,j]
            for i in nodes, j in nodes, m in 1:NV
        ) +
        sum(K(k).vl * dist(i,j) * routes[m].w[i,j]
            for i in nodes, j in nodes, m in 1:NV
        ) +
        sum(K(k).fd * routes[m].u[i]
            for i in nodes, m in 1:NV
        ) +
        sum(K(k).fp * routes[m].z[j]
            for j in nodes, m in 1:NV
        ) +
        0.0000001 * sum(routes[m].ind for m in 1:NV)
    )

    optimize!(res)

    return routes
end

function disaggregate(deli,t)
    results = Dict()

    to_disaggr = []
    for k in K()
        if value.(sum(deli[(k,t)].z[j] for j in K(k).cover)) > 0
            push!(to_disaggr,k)
        end
    end #list vehicles used at period

    for k_dis in to_disaggr
        println("#==KENDARAAN $k_dis, capacity: $(K(k_dis).Q)==#")
        println()

        seed = Vector{Int64}()
        for i in K(k_dis).cover
            if value.(deli[(k_dis,t)].z[i]) > 0
                push!(seed,i)
            end
        end #extract all starting point for vehicle group

        for j in seed
            result = generate(deli,j,k_dis,t)

            results[(j = j,k = k_dis,t = t)] = result[1]

            println()
            println("FINAL RESULT:")
            println("banyak rute: $(length(result[1]))")
            println("kendaraan tidak tergunakan: $(result[3])")
            println()
        end #disggregate to its starting point
    end

    return results
end

function costaggr(R,k_aggr,t_aggr)
    val = (
        sum(
            sum(2 * K(k).vx * dist(i,j) * value(R[(k,t)].p[i,j])
                for i in K(k).cover, j in K(k).cover
            )
            for k in k_aggr, t in t_aggr
        ) +
        sum(
            sum(K(k).vl * dist(i,j) * value(R[(k,t)].o[i,j])
                for i in K(k).cover, j in K(k).cover
            )
            for k in k_aggr, t in t_aggr
        ) +
        sum(
            sum(K(k).fd * value(R[(k,t)].u[i])
                for i in K(k).cover
            )
            for k in k_aggr, t in t_aggr
        ) +
        sum(
            sum(K(k).fp * value(R[(k,t)].z[j])
                for j in K(k).cover
            )
            for k in k_aggr, t in t_aggr
        )
    )

    return val
end

function costdis(routes::Dict,k_dis,t_dis)
    to_cost = filter(p -> first(p).k in k_dis && first(p).t in t_dis,routes)

    val = 0
    for r in to_cost
        for m in last(r)
            #println(first(r))
            #println(first(m))
            cost = sum(K(first(r).k).vx * dist(i,j) * last(m).x[i,j]
                for i in last(m).x.axes[1], j in last(m).x.axes[2]
            ) +
            sum(K(first(r).k).vl * dist(i,j) * last(m).w[i,j]
                for i in last(m).w.axes[1], j in last(m).w.axes[2]
            ) +
            sum(K(first(r).k).fd * last(m).u[i]
                for i in last(m).u.axes[1]
            ) +
            sum(K(first(r).k).fp * last(m).z[j]
                for j in last(m).z.axes[1]
            )

            val += cost
        end
    end

    return val
end
