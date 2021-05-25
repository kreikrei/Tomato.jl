# =========================================================================
#    BASIC DATA GENERATION
# =========================================================================

const vertex_data = Ref{Any}(nothing)
V() = sort!(collect(keys(vertex_data[])))
V(i) = vertex_data[][i]

const vehicle_data = Ref{Any}(nothing)
K() = sort!(collect(keys(vehicle_data[])))
K(k) = vehicle_data[][k]

const period_data = Ref{Any}(nothing)
T() = period_data[]
T(t) = period_data[][t]

const demand_data = Ref{Any}(nothing)
d() = demand_data[]
d(i,t) = demand_data[][i,t]

const distance_data = Ref{Any}(nothing)
dist() = distance_data[]
dist(i,j) = distance_data[][i,j]

function extract!(path::String;f::String) #extract from excel
    xf = XLSX.readxlsx(path) #READ WORKSHEET
    data = Dict{Symbol,DataFrame}() #DATAFRAME DICT

    for sheets in XLSX.sheetnames(xf) #TURN SHEETS INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...) #TRANSFORM TABLE INTO DATAFRAME
        data[Symbol(sheets)] = df #DEFINE THE NAME FROM THE WORKSHEET
    end

    V = process_vertex(data[:vertices])
    K = process_vehicle(data[:vehicles])

    T = collect( #range from starting month for duration
        range(
            last(data[:periods].start),
            length = last(data[:periods].T),
            step = 1
        )
    )

    dist = JuMP.Containers.DenseAxisArray{Float64}(undef, keys(V), keys(V))
    dist .= 999999999
    for i in keys(V), j in keys(V)
        if i != j
            if f == "haversine"
                dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],6378.137)
            else #if f == "euclidean"
                dist[i,j] = euclidean([V[i].x,V[i].y],[V[j].x,V[j].y])
            end
        end
    end

    #=d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:non_negative_demands][:,string.(T)]), #dataset
        Array{String}(data[:non_negative_demands].point), #dims 1
        T #dims 2
    )=#

    d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:demands][:,string.(T)]), #dataset
        Array{String}(data[:demands].point), #dims 1
        T #dims 2
    )

    vertex_data[] = V
    vehicle_data[] = K
    period_data[] = T
    demand_data[] = d
    distance_data[] = dist

    return V,K,T,d,dist
end

function process_vertex(df::DataFrame)
    V = Dict{String,vtx}() #INITIATE VERTICES
    for v in eachrow(df) #ITERATE OVER DATA
        V[v.name] = vtx(
            v.x, v.y, v.MAX, v.MIN, v.START, v.h
        )
    end

    return V
end

function process_vehicle(df::DataFrame)
    K = Dict{String,veh}() #INITIATE VEHICLES
    idx_asal = unique(df.Asal)
    idx_tujuan = unique(df.Tujuan)
    list_moda = unique(df.Moda)

    for m in list_moda
        trayek = filter(p -> p.Moda == m, df)
        unique!(trayek) #trayek yang unik

        Q = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        f = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        g = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        lim = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)

        Q .= 0
        f .= 0
        g .= 0
        lim .= 0

        for r in eachrow(trayek)
            Q[r.Asal,r.Tujuan] = round(r.kapasitas)
            f[r.Asal,r.Tujuan] = round(r.peti_cost)
            g[r.Asal,r.Tujuan] = round(r.trip_cost)
            lim[r.Asal,r.Tujuan] = round(r.limit)
        end

        K[m] = veh(Q,f,g,lim)
    end

    return K
end

const default_optimizer = Ref{Any}(nothing)
set_optimizer!(O) = default_optimizer[] = O
get_optimizer() = default_optimizer[]
reset_optimizer() = default_optimizer[] = nothing
