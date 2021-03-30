# =========================================================================
#    ENVIRONMENT SETTINGS
# =========================================================================

#OPTIMIZER SETTING
const default_optimizer = Ref{Any}(nothing)
set_optimizer!(O) = default_optimizer[] = O
get_optimizer() = default_optimizer[]
reset_optimizer() = default_optimizer[] = nothing
