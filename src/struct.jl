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
    load::Float64 #jumlah peti
    trip::Float64 #trip yang dilakukan
    k::String #jenis kendaraan
    t::Int64 #period

    function delivery(asal,tujuan,load,trip,k,t)
        if !(asal in V()) || !(tujuan in V())
            error("not in area")
            return nothing
        elseif load > K(k).Q[asal,tujuan] * ceil(trip)
            error("overload")
            return nothing
        elseif trip > K(k).lim[asal,tujuan]
            error("overtrip")
            return nothing
        else
            return new(asal,tujuan,load,trip,k,t)
        end
    end
end
