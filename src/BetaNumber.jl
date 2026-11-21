module BetaNumber

using JuMP
using MosekTools
using Dualization

export beta_number, delta_number

mutable struct beta_data
    tbasis
    basis
    moment
end


include("utils.jl")
include("sdp.jl")

end
