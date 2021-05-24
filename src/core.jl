# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================
function original()
    feas = Model(get_optimizer())
    #set_silent(feas)

    set_optimizer_attribute(feas, "TimeLimit", 40)

    @variable(feas, o[i = V(), j = V(), k = K(), t = T()] >= 0, Int)
    @variable(feas, p[i = V(), j = V(), k = K(), t = T()] >= 0, Int)
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

    ori_sol = Vector{delivery}()
    for i in V(), j in V(), k in K(), t in T()
        if value(p[i,j,k,t]) > 0
            new_delivery = delivery(i,j,value(o[i,j,k,t]),value(p[i,j,k,t]),k,t)
            push!(ori_sol, new_delivery)
        end
    end #COLLECT SOLUTION

    return ori_sol
end #ORIGINAL FORMULATION OF IRP

function relaxed()
    feas = Model(get_optimizer())
    set_silent(feas)

    @variable(feas, o[i = V(), j = V(), k = K(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, p[i = V(), j = V(), k = K(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())], Int) #RELAXED INVENTORY

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
    end #COLLECT SOLUTION

    return relaxed_sol
end #SAME AS ORI BUT NO INTEGER DEFINITION

function cluster()
    Q_cluster = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V())
    f_cluster = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V())
    g_cluster = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V())
    lim_cluster = JuMP.Containers.DenseAxisArray{Float64}(undef,V(),V())

    Q_cluster .= 0
    f_cluster .= 0
    g_cluster .= 0
    lim_cluster .= 0

    for i in V(), j in V()
        active_mode = Vector{String}()
        for k in K()
            if K(k).lim[i,j] > 0
                push!(active_mode,k)
            end
        end #CHECK ALL ACTIVE MODE OF TRANSPORT

        if length(active_mode) >= 1
            weight_vec = Dict{String,Float64}()
            total_weight = sum(K(m).lim[i,j] for m in active_mode)
            for m in active_mode
                weight_vec[m] = K(m).lim[i,j]/total_weight
            end #vektor bobot proporsi terhadap total biaya peti

            Q_cluster[i,j] = ceil(sum(weight_vec[m] * K(m).Q[i,j] for m in active_mode))
            f_cluster[i,j] = sum(weight_vec[m] * K(m).f[i,j] for m in active_mode)
            g_cluster[i,j] = sum(weight_vec[m] * K(m).g[i,j] for m in active_mode)
            lim_cluster[i,j] = ceil(sum(weight_vec[m] * K(m).lim[i,j] for m in active_mode))
        end #PROSES SEKIRANYA ADA KENDARAAN DI SEGMEN
    end

    return veh(Q_cluster,f_cluster,g_cluster,lim_cluster)
end #OUTPUTS NEW AGGREGATED VEHICLE

function reduced(cluster::veh)
    feas = Model(get_optimizer())
    set_silent(feas)

    @variable(feas, o[i = V(), j = V(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, p[i = V(), j = V(), t = T()] >= 0) #RELAXED DELIVERY
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())], Int) #RELAXED INVENTORY

    @constraint(feas, [i = V(), t = T()],
        I[i,t-1] - I[i,t] - d(i,t) +
        sum(o[j,i,t] for j in V()) - sum(o[i,j,t] for j in V()) == 0
    ) #Inventory Balance

    @constraint(feas, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    ) #Inventory Max Min

    @constraint(feas, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    ) #Inventory Start

    @constraint(feas, [i = V(), j = V(), t = T()],
        o[i,j,t] <= cluster.Q[i,j] * p[i,j,t]
    ) #quantity to trips proportion from ij

    @constraint(feas, [i = V(), j = V(), t = T()],
        p[i,j,t] <= cluster.lim[i,j]
    ) #limit of trips from ij

    @objective(feas, Min,
        sum(
            cluster.f[i,j] * o[i,j,t] + cluster.g[i,j] * p[i,j,t]
            for i in V(), j in V(), t in T()
        ) +
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        )
    )

    optimize!(feas)

    deliveries = Vector{delivery}()
    for i in V(), j in V(), t in T()
        if value(p[i,j,t]) > 0
            dis = Model(get_optimizer())
            set_silent(dis)

            @variable(dis, load[k = K()] >= 0, Int)
            @variable(dis, trip[k = K()] >= 0, Int)

            @constraint(dis, sum(load) >= value(o[i,j,t]))
            @constraint(dis, [k = K()], load[k] <= K(k).Q[i,j] * trip[k])
            @constraint(dis, [k = K()], trip[k] <= K(k).lim[i,j])

            @objective(dis, Min,
                sum(K(k).f[i,j] * load[k] + K(k).g[i,j] * trip[k] for k in K())
            ) #sum of all transport mode

            optimize!(dis)

            for k in K()
                if value(trip[k]) > 0
                    new_delivery = delivery(i,j,value(load[k]),value(trip[k]),k,t)
                    push!(deliveries,new_delivery)
                end
            end
        end
    end #DISAGGREGATE SOLUTION FOR EACH SEGMENT

    return deliveries
end #COMPUTE AGGREGATED SOLUTION

function roundUp(relaxed_sol::Vector{delivery})
    deliveries = Vector{delivery}()

    for r in relaxed_sol
        new_delivery = delivery(r.asal,r.tujuan,r.load,ceil(r.trip),r.k,r.t)
        push!(deliveries, new_delivery)
    end

    return deliveries
end #ALTERNATIVE FOR FEASIBLE SOL

function splatBuild(relaxed_sol::Vector{delivery})
    deliveries = Vector{delivery}()
    unique_asal = unique([p.asal for p in relaxed_sol])
    unique_tujuan = unique([p.tujuan for p in relaxed_sol])
    unique_t = unique([p.t for p in relaxed_sol])

    for i in unique_asal, j in unique_tujuan, t in unique_t
        active_delivery = filter(p -> p.asal == i && p.tujuan == j && p.t == t,relaxed_sol)
        aggregate = sum([r.load for r in active_delivery]) #collect aggregate delivery

        if aggregate > 0
            dis = Model(get_optimizer())
            set_silent(dis)

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

function optimizeVolume(deliveries::Vector{delivery})
    OV = Model(get_optimizer())
    set_silent(OV)

    @variable(OV, o[r = deliveries] >= 0)
    @variable(OV, I[i = V(), t = vcat(first(T()) - 1, T())], Int)

    @constraint(OV, [i = V(), t = T()],
        I[i,t-1] - I[i,t] - d(i,t) +
        sum(o[r] for r in filter(p -> p.tujuan == i && p.t == t,deliveries)) -
        sum(o[r] for r in filter(p -> p.asal == i && p.t == t,deliveries)) == 0
    ) #Inventory Balance

    @constraint(OV, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    ) #Inventory Max Min

    @constraint(OV, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    ) #Inventory Start

    @constraint(OV, [r = deliveries],
        o[r] <= K(r.k).Q[r.asal,r.tujuan] * r.trip
    ) #quantity to trips proportion from ij

    @objective(OV, Min,
        sum(
            o[r] * K(r.k).f[r.asal,r.tujuan] + r.trip * K(r.k).g[r.asal,r.tujuan]
            for r in deliveries
        ) +
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        )
    )

    optimize!(OV)

    replacement = Vector{delivery}()
    for r in deliveries
        new_load = value(o[r])
        new_trip = ceil(new_load/K(r.k).Q[r.asal,r.tujuan]) #revise trip

        if new_trip > 0
            new_delivery = delivery(r.asal, r.tujuan, new_load, new_trip, r.k, r.t)
            push!(replacement,new_delivery)
        end
    end #REVISE TRIP NUMBER

    return replacement
end #reoptimize the volume from trips

function iterOV(deliveries::Vector{delivery})
    mem = length(deliveries)
    deliveries = optimizeVolume(deliveries)

    while mem != length(deliveries)
        mem = length(deliveries)
        deliveries = optimizeVolume(deliveries)
    end

    return deliveries
end #iterate optimize volume until no improvement

function deliveryCost(deliveries)
    if isa(deliveries,Vector{delivery})
        total_cost = 0

        for r in deliveries
            indv_cost = r.load * K(r.k).f[r.asal,r.tujuan] + r.trip * K(r.k).g[r.asal,r.tujuan]
            total_cost += indv_cost
        end

        return total_cost
    else
        indv_cost = deliveries.load * K(deliveries.k).f[deliveries.asal,deliveries.tujuan] + deliveries.trip * K(deliveries.k).g[deliveries.asal,deliveries.tujuan]

        return indv_cost
    end
end #calculate cost of deliveries

function inventoryCost(deliveries::Vector{delivery})
    total_cost = 0
    I = stockLevels(deliveries)
    total_cost = sum(V(i).h * I[i,t] for i in V(), t in T())

    return total_cost
end #calculate inventory costs implied by deliveries

totalCost(deliveries) = deliveryCost(deliveries) + inventoryCost(deliveries)

function stockLevels(deliveries::Vector{delivery})
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
end #compute stock levels

function feasibility(deliveries::Vector{delivery})
    I = stockLevels(deliveries)
    infeasibilities = Vector{NamedTuple}()

    for i in V(), t in T()
        if I[i,t] > V(i).MAX || I[i,t] < V(i).MIN
            push!(infeasibilities,(i=i,t=t))
        end
    end

    if !isempty(infeasibilities)
        println("infeasible")
    end

    return infeasibilities
end #track infeasibilities in inventory level

function destroy(deliveries::Vector{delivery})
    c = rand(V()) #random combination of random length

    destroyed = deepcopy(deliveries) #make independent duplicate
    removed = filter(p -> p.tujuan == c || p.asal == c, destroyed) #connected to c set
    for r in removed
        deleteat!(destroyed, findall(x -> x == r, destroyed))
    end #remove all

    return destroyed,c,removed
end #make the solution infeasible

function repair(destroyed::Vector{delivery},c::String,removed::Vector{delivery})
    feas = Model(get_optimizer())
    set_silent(feas)

    o★ = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V(),K(),T())
    p★ = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V(),K(),T())

    o★ .= 0
    p★ .= 0

    for i in V(), j in V(), k in K(), t in T()
        existing = filter(p -> p.asal==i && p.tujuan==j && p.k==k && p.t==t, destroyed)
        if !isempty(existing)
            o★[i,j,k,t] = round(existing[1].load)
            p★[i,j,k,t] = round(existing[1].trip)
        end
    end #transform existing deliveries into matrix

    to_repair = filter(p ->  p.asal == c || p.tujuan == c,all_edges())
    #println(to_repair)

    @variable(feas, o[e = to_repair, k = K(), t = T()] >= 0, Int)
    @variable(feas, p[e = to_repair, k = K(), t = T()] >= 0, Int)
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())], Int)

    @constraint(feas, balance[i = V(), t = T()],
        I[i,t-1] - I[i,t] - d(i,t) + sum(
            sum(o★[j,i,k,t] for j in V()) -
            sum(o★[i,j,k,t] for j in V())
            for k in K()
        ) + sum(
            sum(o[e,k,t] for e in to_repair if e.tujuan == i) -
            sum(o[e,k,t] for e in to_repair if e.asal == i)
            for k in K()
        ) == 0
    ) #inventory balance

    @constraint(feas, [i = V(), t = T()],
        V(i).MIN <= I[i,t] <= V(i).MAX
    ) #Inventory Max Min

    @constraint(feas, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    ) #Inventory Start

    @constraint(feas, [e = to_repair, k = K(), t = T()],
        o[e,k,t] <= K(k).Q[e.asal,e.tujuan] * p[e,k,t]
    ) #load to trips proportion from ij

    @constraint(feas, [e = to_repair, k = K(), t = T()],
        p[e,k,t] <= K(k).lim[e.asal,e.tujuan]
    ) #limit of trips from ij

    @objective(feas, Min,
        sum(
            K(k).g[e.asal,e.tujuan] * p[e,k,t] + K(k).f[e.asal,e.tujuan] * o[e,k,t]
            for e in to_repair, k in K(), t in T()
        ) +
        sum(
            V(i).h * I[i,t]
            for i in V(), t in T()
        )
    )

    #=@objective(feas, Min,
        sum(
            K(k).g[e.asal,e.tujuan] * p[e,k,t]
            for e in to_repair, k in K(), t in T()
        )
    )=#

    set_optimizer_attribute(feas, "Presolve", 2)
    set_optimizer_attribute(feas, "TimeLimit", 0.1)
    set_optimizer_attribute(feas, "MIPFocus", 1)

    optimize!(feas)

    new_deliveries = Vector{delivery}()
    for e in to_repair, k in K(), t in T()
        if value(p[e,k,t]) > 0
            new_delivery = delivery(
                e.asal,e.tujuan,
                round(value(o[e,k,t])),
                ceil(value(o[e,k,t])/K(k).Q[e.asal,e.tujuan]),
                k,t
            )
            push!(new_deliveries,new_delivery)
        end
    end #EXTRACT ADDITIONAL SOLUTION

    append!(destroyed,new_deliveries)

    return destroyed
end #make the solution feasible, if not optimal

function search(x::Vector{delivery})
    x_best = deepcopy(x)
    counter = 0

    limit = 200
    while counter < limit && totalCost(x_best) > 3.5e10
        destroyed,c,removed = destroy(x)
        x_temp = repair(destroyed,c,removed)

        if isempty(feasibility(x_temp))
            if totalCost(x_temp) < totalCost(x)
                x = deepcopy(x_temp)
            end

            if totalCost(x_temp) < totalCost(x_best)
                x_best = deepcopy(x_temp)
                println("new best cost: $(totalCost(x_best)) after $counter iter")
                counter = 0 #resest counter
            else
                counter += 1 #no improvement
                println("counter: $counter")
            end
        end
    end

    return x_best
end

function stockouts(deliveries)
    I = stockLevels(deliveries)

    for i in V(), t in T()
        if I[i,t] >= V(i).MIN
            I[i,t] = 0
        else
            I[i,t] = V(i).MIN - I[i,t]
        end
    end

    return I
end

function overstocks(deliveries)
    I = stockLevels(deliveries)

    for i in V(), t in T()
        if I[i,t] <= V(i).MAX
            I[i,t] = 0
        else
            I[i,t] = I[i,t] - V(i).MAX
        end
    end

    return I
end

function available_supply(deliveries)
    I = stockLevels(deliveries)

    for i in V(), t in T()
        if I[i,t] > V(i).MIN
            I[i,t] = I[i,t] - V(i).MIN
        else
            I[i,t] = 0
        end
    end

    return I
end

function available_vaccancy(deliveries)
    I = stockLevels(deliveries)

    for i in V(), t in T()
        if I[i,t] < V(i).MAX
            I[i,t] = V(i).MAX - I[i,t]
        else
            I[i,t] = 0
        end
    end

    return I
end
