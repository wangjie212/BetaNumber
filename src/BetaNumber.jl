module BetaNumber

using JuMP
using MosekTools

export beta_number

include("utils.jl")
include("sdp.jl")

end
