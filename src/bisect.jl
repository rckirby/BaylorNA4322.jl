using Printf

"""
    bisect(f, a, b, eps, maxit, output=false)

    Uses the bisection method to find an `x` between `a` and `b` such that 
    the univariate function `f` is close to zero at `x`.  Stopping criteria
    is that the width of the interval obtained by bisection is less than 
    `eps`.  Takes no more than `maxit` steps.

    If `output` is set to `true`, then diagnostic information is displayed
    at each iteration.

    returns a triple of values (converged, its, c) indicating
    whether or not the method converged (a boolean), the number of iterations
    taken, and the final approximation obtained to the root.

    This will raise an exception if at any point it cannot determine that
    the interval it is working on contains a root (i.e. if `f` has the same 
    sign at the endpoints.
"""
function bisect(f, a, b, eps, maxit, output=false)
    fa = f(a)
    fb = f(b)
    for (x, fx) = [(a,fa),(b,fb)]
        if fx == 0.0
            return true, x, 0
        end
    end
    if fa * fb > 0
        throw(ErrorException("root not bracketed"))
    end
    
    i = 0
    
    while i < maxit
        @assert(fb*fa < 0)
        if output
            @printf "%5d\t%0.6e\t%0.6e\t%0.6e\t%0.6e\n" i a fa b fb
        end
        c = 0.5 * (a + b)
        fc = f(c)
        i += 1
        if fa * fc < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
        if b-a < eps
            return true, i, c
        end
    end
    return false, i, c
end
