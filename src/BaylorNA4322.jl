module BaylorNA4322

greet() = print("Hello World!")

include("float.jl")
include("roots.jl")
include("interp.jl")
export FloatingPointNumbers
export RootFinding
export Interpolation
end # module
