module BetaNumber

using JuMP
using MosekTools

export beta_number

mutable struct beta_data
    tbasis
    basis
    moment
end


include("utils.jl")
include("sdp.jl")

end
