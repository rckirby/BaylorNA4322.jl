module Interpolation
include("neville.jl")
include("newtlag.jl")
export neville
export newton_divided_diff_tableau, newton_divided_diff, newton_eval
end
