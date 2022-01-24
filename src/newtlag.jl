using LinearAlgebra  # for "diag"

"""
    newton_divided_diff_tableau(xi, fi)

    Computes the divided difference tableau using nodes xi
    and function values fi.  Returns a dense matrix of values.

    The coefficients of the Newton form of the polynomial will
    lie on the diagonal of this matrix.
"""
function newton_divided_diff_tableau(xi, fi)
    N = length(xi)
    tableau = zeros(N, N)
    tableau[:, 1] = fi[:]
    for j=2:N
        for i=j:N
            tableau[i, j] = (tableau[i, j-1] - tableau[i-1, j-1])/(xi[i]-xi[i-j+1])
        end
    end
    
    return tableau
end


"""
    newton_divided_diff(xi, fi)

    Computes the Newton coefficients of an interpolating polynomial
    using divided differences.  Only returns the coefficients, does not
    store the entire tableau.
"""
function newton_divided_diff(xi, fi)
    # *just* computes the Newton form coefficients of the Lagrange interpolating
    # polynomial
    N = length(xi)
    c = copy(fi)
    for j=2:N
        for i=N:-1:j
            c[i] = (c[i] - c[i-1])/(xi[i]-xi[i-j+1])
        end
    end
    return c
end


"""
    newton_eval(xi, c, x)

    Evaluates the Newton form of an interpolating polynomial at a point x.
    The Lagrange nodes are stored in xi, and the Newton coefficients in c.
"""
function newton_eval(xi, c, x)
    N = length(c)
    val = c[N]
    for i=N-1:-1:1
        val = val*(x-xi[i]) + c[i]
    end
    return val
end

"""Creates a vector that repeats each entry"""
function twice(xi)
    n = length(xi)
    xx = zeros(2*n)
    for i=1:n
        xx[2*i-1] = xi[i]
        xx[2*i]  = xi[i]
    end
    return xx
end

"""
    newton_hermite_divided_diff(xi, fi, fip)

    Computes the divided difference tableau using nodes xi,
    function values fi, and derivatives fip.  

    Returns the vector of repeated nodes and the
    Newton coefficients of the polynomial.
"""
function newton_hermite_divided_diff(xi, fi, fip)
    n = length(xi)
    @assert n == length(fi)
    @assert n == length(fip)
    N = 2*n

    tableau = zeros(N, N)
    # First column has repeated function values
    for i=1:n
        tableau[2*i-1, 1] = fi[i]
        tableau[2*i, 1]= fi[i]
    end

    # Second column alternates derivatives and differences
    tableau[2, 2] = fip[1]
    for i=2:n
        tableau[2*i-1, 2] = (fi[i] - fi[i-1]) / (xi[i] - xi[i-1])
        tableau[2*i, 2] = fip[i]
    end

    xx = twice(xi)
    # Now, differences as "normal"!
    for j=3:N
        for i=j:N
            tableau[i, j] = (tableau[i, j-1] - tableau[i-1, j-1])/(xx[i]-xx[i-j+1])
        end
    end
    
    return xx, diag(tableau)
end


using AutoGrad
f(x) = exp(-2x)
fp = grad(f)
xs = [0., 1.]
fs = f.(xs)
fps = fp.(xs)

xx, c = newton_hermite_divided_diff(xs, fs, fps)

print(newton_eval(xx, c, 0.5))

