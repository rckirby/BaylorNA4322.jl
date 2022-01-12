using Printf


"""
    newton(f, fp, x0, eps, maxit, output=false)

    Given a function f and its derivative fp, performs Newton iteratin
    starting from x0 to approximate a root of f(x).  The function terminates
    when |f(x)| < eps or maxit iterations have been taken.

    returns a triple of values (converged, its, c) indicating
    whether or not the method converged (a boolean), the number of iterations
    taken, and the final approximation obtained to the root.
"""
function newton(f, fp, x0, eps, maxit, output=false)
    x = x0
    i = 0
    fx = f(x)
    if output
        @printf "  Its   x               f(x)            \n"
        @printf "==========================================\n"
    end
    while i < maxit
        if abs(fx) < eps
            return true, i, x
        end
        fx = f(x)
        fpx = fp(x)
        x -= fx / fpx
        i += 1
        if output
            @printf "%5d\t%0.8e\t%0.8e\n" i x fx
        end
    end
    return false, i, x
end
