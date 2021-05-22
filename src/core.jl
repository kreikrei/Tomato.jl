# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================
function cluster()
    Q_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    f_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    g_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    lim_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())

    Q_cluster .= 0
    f_cluster .= 0
    g_cluster .= 0
    lim_cluster .= 0

    for i in V(), j in V()
        active_mode = Vector{String}()
        for k in K()
            if K(k).Q[i,j] > 0
                push!(active_mode,k)
            end
        end #CHECK ALL ACTIVE MODE OF TRANSPORT

        if length(active_mode) >= 1
            weight_vec = Dict{String,Float64}()
            total_weight = sum(K(m).f[i,j] for m in active_mode)
            for m in active_mode
                weight_vec[m] = (K(m).f[i,j])/total_weight
            end #vektor bobot proporsi terhadap total biaya peti

            Q_cluster[i,j] = round(sum(weight_vec[m] * K(m).Q[i,j] for m in active_mode))
            f_cluster[i,j] = round(sum(weight_vec[m] * K(m).f[i,j] for m in active_mode))
            g_cluster[i,j] = round(sum(weight_vec[m] * K(m).g[i,j] for m in active_mode))
            lim_cluster[i,j] = ceil(sum(weight_vec[m] * K(m).lim[i,j] for m in active_mode))
        end #PROSES SEKIRANYA ADA KENDARAAN DI SEGMEN
    end

    return Q_cluster,f_cluster,g_cluster,lim_cluster
end

function reduced(Q,f,g,lim)
    feas = Model(get_optimizer())

    @variable(feas, o[i = V(), j = V(), t = T()] >= 0)
    @variable(feas, p[i = V(), j = V(), t = T()] >= 0)
    @variable(feas, S[i = V(), t = T()] >= 0, Int)
    @variable(feas, I[i = V(), t = vcat(first(T()) - 1, T())], Int)

    @constraint(feas, [i = V(), t = T()],
        I[i,t-1] + S[i,t] + sum(o[j,i,t] for j in V()) - sum(o[i,j,t] for j in V()) ==
        d(i,t) + I[i,t]
    ) #inventory balance

    @constraint(feas, [i = V(), t = T()],
        V(i).MIN <= I[i,t]
    ) #inventory max min

    @constraint(feas, [i = V(), t = T()],
        I[i,t] <= V(i).MAX
    ) #inventory max min

    @constraint(feas, [i = V()],
        I[i,first(T()) - 1] == V(i).START
    ) #inventory start

    @constraint(feas, [i = V(), j = V(), t = T()],
        o[i,j,t] <= Q[i,j] * p[i,j,t]
    ) #segment usage limit

    @constraint(feas,[i = V(), j = V(), t = T()],
        p[i,j,t] <= lim[i,j]
    )

    @constraint(feas, [i = V(), t = T()],
        S[i,t] <= 999999999999 * V(i).prod
    )

    @objective(feas, Min,
        sum(
            f[i,j] * o[i,j,t] + g[i,j] * p[i,j,t]
            for i in V(), j in V(), t in T()
        ) +
        sum(
            500000000 * S[i,t] + V(i).h * I[i,t]
            for i in V(), t in T()
        )
    )

    optimize!(feas)

    sol = Vector{aggr}()
    for i in V(), j in V(), t in T()
        if value(o[i,j,t]) > 0
            push!(sol, aggr(i,j,round(value(o[i,j,t]))))
        end
    end

    return feas,sol
end

function dissect()

end
