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
