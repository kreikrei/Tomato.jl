# =========================================================================
#    BASIC STRUCTURES
# =========================================================================

struct vtx
    x::Float64 #xcoor
    y::Float64 #ycoor
    MAX::Int64 #max inventory level
    MIN::Int64 #min inventory level
    START::Int64 #starting inventory level
    h::Float64 #inv cost per unit
end

struct veh
    Q::JuMP.Containers.DenseAxisArray #kapasitas segmen
    f::JuMP.Containers.DenseAxisArray #peti_cost
    g::JuMP.Containers.DenseAxisArray #trip_cost
    lim::JuMP.Containers.DenseAxisArray #usage_limit
end

struct delivery
    asal::String
    tujuan::String
    q::Int64 #tidak peduli berapa pengiriman
    k::String #jenis kendaraan
    t::Int64 #period
end

struct disaggr
    asal::String
    tujuan::String
    kendaraan::String
    q::Int64
end
