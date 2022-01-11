module FloatingPointNumbers

export make_fl_float

struct FloatingPoint{N}
    mantissa::Int
    characteristic::Int
end

"""
    make_FP(x, N)

    makes an N-digit `FloatingPoint{N}` number from some Julia floating
    point type.  Returns a number of the form +/- M * 10^n,
    where the mantissa satisfies 0.1 <= M < 1.

    Future versions may account for overflow/underflow
"""
function make_FP(x, N)
    # get the number in the form x = +/- M * 10^n,
    # where 0.1 <= M < 1
    
    # check for underflow, set number to 0
    if abs(x) < 1.e-12
        return FloatingPoint{N}(0, 0)
    end
    n = Int(floor(log10(abs(x)) + 1))
    absx = abs(x)
    sgnx = sign(x)
    mant = Int(floor(absx*10.0^(N-n)))
    return FloatingPoint{N}(sgnx*mant, n)
end

Base.show(io::IO, x::FloatingPoint) = print(io, "Fl(", x.mantissa>0 ? "" : "-", "0.", abs(x.mantissa), " x 10^", x.characteristic, ")")

"Converts a FloatingPoint{N} number back to a Julia float"
FP2float(x, N) = x.mantissa * 10.0^(x.characteristic-N)

"Creates a pair of functions for converting to/from the FloatingPoing{N} representation."
make_fl_float(N) = (x->make_FP(x, N), x->FP2float(x, N))


# overload arithmetic to DTRT with FloatingPoint{N} numbers.
import Base: +, -, *, /, ^, abs, sqrt, <,<=,>,>=,==,!=
for op in (:+, :-, :*, :/, :^)
   @eval begin
        $op(x::FloatingPoint{N}, y::FloatingPoint{N}) where N = make_FP(($op)(FP2float(x,N),FP2float(y,N)),N)
   end
end

# Lots of other unary operations (trig functions, exp, etc) could
# be included here but aren't.
(-)(x::FloatingPoint{N}) where N = FloatingPoint{N}(-x.mantissa, x.characteristic)
sqrt(x::FloatingPoint{N}) where N = make_FP(sqrt(FP2float(x,N)),N)
abs(x::FloatingPoint{N}) where N = FloatingPoint{N}(abs(x.mantissa), x.characteristic)

for op in [:<,:<=,:>,:>=,:(==),:!=]
    @eval begin
        $op(x::FloatingPoint{N}, y::FloatingPoint{N}) where N = ($op)(FP2float(x,N),FP2float(y,N))
    end
end

end
