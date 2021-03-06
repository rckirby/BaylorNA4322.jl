module Interpolation
include("neville.jl")
include("newtlag.jl")
include("spline.jl")
export neville, neville_table
export newton_divided_diff_tableau, newton_divided_diff, newton_eval, newton_hermite_divided_diff
export compute_natural_spline, evaluate_spline, evaluate_spline_derivative
end
