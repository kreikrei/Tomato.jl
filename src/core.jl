# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================
function relaxed()
    feas = Model(get_optimizer())

    @variable(feas, o[i = V(), j = V(), k = K(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, p[i = V(), j = V(), k = K(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())], Int)

    @constraint(feas, [i = V(), t = T()],
        I[i,t-1] - I[i,t] - d(i,t) +
        sum(sum(o[j,i,k,t] for j in V()) - sum(o[i,j,k,t] for j in V()) for k in K()) == 0
    ) #Inventory Balance

    @constraint(feas, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    ) #Inventory Max Min

    @constraint(feas, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    ) #Inventory Start

    @constraint(feas, [i = V(), j = V(), k = K(), t = T()],
        o[i,j,k,t] <= K(k).Q[i,j] * p[i,j,k,t]
    ) #quantity to trips proportion from ij

    @constraint(feas, [i = V(), j = V(), k = K(), t = T()],
        p[i,j,k,t] <= K(k).lim[i,j]
    ) #limit of trips from ij

    @objective(feas, Min,
        sum(
            K(k).f[i,j] * o[i,j,k,t] + K(k).g[i,j] * p[i,j,k,t]
            for i in V(), j in V(), k in K(), t in T()
        ) +
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        )
    )

    optimize!(feas)

    relaxed_sol = Vector{delivery}()
    for i in V(), j in V(), k in K(), t in T()
        if value(p[i,j,k,t]) > 0
            new_delivery = delivery(i,j,value(o[i,j,k,t]),value(p[i,j,k,t]),k,t)
            push!(relaxed_sol, new_delivery)
        end
    end

    return relaxed_sol
end

function roundUp(relaxed_sol)
    deliveries = Vector{delivery}()

    for r in relaxed_sol
        println(r)
        new_delivery = delivery(r.asal,r.tujuan,r.load,ceil(r.trip),r.k,r.t)
        push!(deliveries, new_delivery)
    end

    return deliveries
end #ALTERNATIVE FOR FEASIBLE SOL

function splatBuild(relaxed_sol)
    deliveries = Vector{delivery}()
    unique_asal = unique([p.asal for p in relaxed_sol])
    unique_tujuan = unique([p.tujuan for p in relaxed_sol])
    unique_t = unique([p.t for p in relaxed_sol])

    for i in unique_asal, j in unique_tujuan, t in unique_t
        active_delivery = filter(p -> p.asal == i && p.tujuan == j && p.t == t,relaxed_sol)
        aggregate = sum([r.load for r in active_delivery]) #collect aggregate delivery

        if aggregate > 0
            dis = Model(get_optimizer())

            @variable(dis, o[k = K()] >= 0, Int)
            @variable(dis, p[k = K()] >= 0, Int)

            @constraint(dis, sum(o) == aggregate)
            @constraint(dis, [k = K()], o[k] <= K(k).Q[i,j] * p[k])
            @constraint(dis, [k = K()], p[k] <= K(k).lim[i,j])

            @objective(dis, Min, sum(K(k).f[i,j] * o[k] + K(k).g[i,j] * p[k] for k in K()))

            optimize!(dis)

            for k in K()
                if value(p[k]) > 0
                    new_delivery = delivery(i,j,value(o[k]),value(p[k]),k,t)
                    push!(deliveries,new_delivery)
                end
            end
        end
    end

    return deliveries
end #ALTERNATIVE FOR FEASIBLE SOL    

function deliveryCost(deliveries)
    total_cost = 0

    for r in deliveries
        indv_cost = r.load * K(r.k).f[r.asal,r.tujuan] + r.trip * K(r.k).g[r.asal,r.tujuan]
        total_cost += indv_cost
    end

    return total_cost
end #calculate cost of deliveries

function inventoryCost(deliveries)
    total_cost = 0
    I = stockLevels(deliveries)
    total_cost = sum(V(i).h * I[i,t] for i in V(), t in T())

    return total_cost
end #calculate inventory costs implied by deliveries

function stockLevels(deliveries)
    I = JuMP.Containers.DenseAxisArray{Int64}(undef, V(), vcat(first(T()) - 1, T()))
    for i in V()
        I[i,first(T()) - 1] = V(i).START
    end #fill starting inventory level

    for i in V(), t in T()
        IN = filter(p -> p.tujuan == i && p.t == t,deliveries)
        OUT = filter(p -> p.asal == i && p.t == t,deliveries)

        I[i,t] = I[i,t-1] + sum([r.load for r in IN]) - sum([r.load for r in OUT]) - d(i,t)
    end

    return I
end
