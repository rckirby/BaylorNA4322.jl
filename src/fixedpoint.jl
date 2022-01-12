using Printf

"""
    fpiter(f, x0, eps, maxit, output=false)

    Performs simple fixed point iteration to search for a point `x` such that
    f(x) is (approximately) x.  The interation is x_{n+1} = f(x_n).

    The algorithm starts from the initial guess x0 and continues
    until subsequent iterates are within `eps` of each other. 

    If `output` is set to `true`, then diagnostic information is displayed
    at each iteration.

    returns a triple of values (converged, its, c) indicating
    whether or not the method converged (a boolean), the number of iterations
    taken, and the final approximation obtained to the root.
"""
function fpiter(f, x0, eps, maxit, output=false)
    x = x0
    i = 1
    if output
        @printf "  Its   x               f(x)            \n"
        @printf "==========================================\n"
    end
    while i < maxit
        fx = f(x)
        if output
            @printf "%d\t%0.6e\t%0.6e\n" i x fx
        end
        # Absolute stopping criteria
        if abs(x-fx) < eps
            return true, i, x
        else
            x = fx
            i += 1
        end            
    end
    return false, i, x
end


"""
    fp_root_finder(f, x0, eps, maxit, output=false)

    Attempts to find a point x where f(x) is approximately zero by finding
    a fixed point of the function g(x) = f(x) + x.

    See `fpiter` for further documentation
"""
function fp_root_finder(f, x0, eps, maxit, output=false)
    return fpiter(x->f(x)+x, x0, eps, maxit, output)
end
