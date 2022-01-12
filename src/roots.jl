module RootFinding

include("bisect.jl")
include("fixedpoint.jl")

export bisect
export fpiter, fp_root_finder
end
