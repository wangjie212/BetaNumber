module BetaNumber

using JuMP
using MosekTools
using COSMO

export beta_number

include("utils.jl")
include("sdp.jl")

end
