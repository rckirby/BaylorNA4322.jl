module Interpolation
include("neville.jl")
include("newtlag.jl")
include("spline.jl")
export neville
export newton_divided_diff_tableau, newton_divided_diff, newton_eval
export compute_natural_spline, evaluate_spline, evaluate_spline_derivative
end
