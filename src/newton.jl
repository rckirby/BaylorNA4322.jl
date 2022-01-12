using Printf

"""
    newton(f, fp, x0, eps, maxit, output=false)

    Given a function f and its derivative fp, performs Newton iteration
    starting from x0 to approximate a root of f(x).  The iteration terminates
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


"""
    secant(f, x0, x1, eps, maxit, output=false)

    Given a function f and two initial guesses x0 and x1, 
    performs secant iteration to approximate a root of f(x).  
    The iteration terminates when |f(x)| < eps or maxit iterations have been 
    taken.

    returns a triple of values (converged, its, c) indicating
    whether or not the method converged (a boolean), the number of iterations
    taken, and the final approximation obtained to the root.
"""
function secant(f, x0, x1, eps, maxit, output=false)
    xx0 = x0
    xx1 = x1
    f0 = f(xx0)
    f1 = f(xx1)
    i = 0
    if output
        @printf "  Its   x               f(x)            \n"
        @printf "==========================================\n"
    end    
    while i < maxit
        if abs(f1) < eps
            return true, i, xx1
        end
        xx2 = (xx0*f1-xx1*f0)/(f1-f0)
        f2 = f(xx2)
        i+=1
        if output
            @printf "%5d\t%0.8e\t%0.8e\n" i xx2 f2
        end
        xx0 = xx1
        xx1 = xx2
        f0 = f1
        f1 = f2
    end
    return false, i, xx1
end


"""
    simplified_newton(f, fp0, x0, eps, maxit, output=false)

    Given a function f and a constant slope fp0, performs simplied Newton/
    constant slopeiteration starting from x0 to approximate a root of f(x).
    The iteration terminates when |f(x)| < eps or maxit iterations have been 
    taken.

    returns a triple of values (converged, its, c) indicating
    whether or not the method converged (a boolean), the number of iterations
    taken, and the final approximation obtained to the root.
"""

function simplified_newton(f, fp0, x0, eps, maxit, output=false)
    # Current guess and iteration count
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
        # evaluate derivative
        fx = f(x)
        # perform simplified Newton update
        x -= fx / fp0
        i += 1
        if output
            @printf "%5d\t%0.8e\t%0.8e\n" i x fx
        end
    end
    return false, i, x
end
