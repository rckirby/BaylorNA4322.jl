module RootFinding

include("bisect.jl")
include("fixedpoint.jl")
include("newton.jl")

export bisect
export fpiter, fp_root_finder
export newton, secant, simplified_newton
end
