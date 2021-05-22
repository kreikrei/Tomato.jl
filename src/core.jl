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

    return feas
end

function roundUp(relaxed_sol)


    return deliveries
end

function splatBuild(relaxed_sol)


    return deliveries
end

function cost(deliveries)


    return total_cost
end
