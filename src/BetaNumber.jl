module BetaNumber

using JuMP
using MosekTools
using Graphs

export beta_number

include("utils.jl")
include("sparsity.jl")
include("sdp.jl")

end
