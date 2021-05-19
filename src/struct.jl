# =========================================================================
#    BASIC STRUCTURES
# =========================================================================

struct vtx
    x::Float64 #xcoor
    y::Float64 #ycoor
    MAX::Int64 #max inventory level
    MIN::Int64 #min inventory level
    START::Int64 #starting inventory level

    #COSTS
    h::Float64 #inv cost per unit
end

struct veh
    #IDENTIFIERS
    id::String

    #PARAM
    cap::JuMP.Containers.DenseAxisArray #kapasitas segmen
    f::JuMP.Containers.DenseAxisArray #peti_cost
    g::JuMP.Containers.DenseAxisArray #dist * jarak_cost
    opsnet::JuMP.Containers.DenseAxisArray #0-1
end

struct aggr
    asal::String
    tujuan::String
    q::Int64 #tidak peduli berapa pengiriman
end

struct disaggr
    asal::String
    tujuan::String
    kendaraan::String
    q::Int64
end
