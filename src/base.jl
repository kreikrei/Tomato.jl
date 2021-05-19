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

    d = JuMP.Containers.DenseAxisArray(
        Array{Float64}(data[:demands][:,string.(T)]), #dataset
        Array{String}(data[:demands].point), #dims 1
        T #dims 2
    )

    vertex_data[] = V
    vehicle_data[] = K
    period_data[] = T
    demand_data[] = d

    return V,K,T,d
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

        kapasitas = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        biayapeti = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        biayajarak = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)
        operasi = JuMP.Containers.DenseAxisArray{Int64}(undef,idx_asal,idx_tujuan)

        kapasitas .= 0
        biayapeti .= 0
        biayajarak .= 0
        operasi .= 0

        for r in eachrow(trayek)
            kapasitas[r.Asal,r.Tujuan] = r.kapasitas
            biayapeti[r.Asal,r.Tujuan] = r.peti_cost
            biayajarak[r.Asal,r.Tujuan] = r.jarak_cost
            operasi[r.Asal,r.Tujuan] = 1
        end

        K[m] = veh(m,kapasitas,biayapeti,biayajarak,operasi)
    end

    return K
end
