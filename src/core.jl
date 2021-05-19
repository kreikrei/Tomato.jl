# =========================================================================
#    CORE FUNCTIONS AND MECHANISMS
# =========================================================================
function cluster()
    Q_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    f_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    g_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())
    ops_cluster = JuMP.Containers.DenseAxisArray{Int64}(undef,V(),V())

    Q_cluster .= 0
    f_cluster .= 0
    g_cluster .= 0
    ops_cluster .= 0

    for i in V(), j in V()
        active_mode = Vector{String}()
        for k in K()
            if K(k).opsnet[i,j] == 1
                push!(active_mode,k)
            end
        end #CHECK ALL ACTIVE MODE OF TRANSPORT

        if length(active_mode) >= 1
            weight_vec = Dict{String,Float64}()
            total_weight = sum(K(m).f[i,j] for m in active_mode)
            for m in active_mode
                weight_vec[m] = K(m).f[i,j]/total_weight
            end #vektor bobot proporsi terhadap total biaya peti

            Q_cluster[i,j] = round(sum(weight_vec[m] * K(m).cap[i,j] for m in active_mode))
            f_cluster[i,j] = round(sum(weight_vec[m] * K(m).f[i,j] for m in active_mode))
            g_cluster[i,j] = round(sum(weight_vec[m] * K(m).g[i,j] for m in active_mode))
            ops_cluster[i,j] = 1 #MASUKIN KE MATRIKS
        end #PROSES SEKIRANYA ADA KENDARAAN DI SEGMEN
    end

    return Q_cluster,f_cluster,g_cluster,ops_cluster
end
