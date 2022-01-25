using LinearAlgebra

"""
    compute_natural_spline(xi, fi)

    xi is an array of n+1 abscissae
    fi is an array of n+1 ordinates

    Constructs a cubic spline interpolant of this data
    using n pieces.

    The spline is represented as a 4 x n array giving the coefficients
    on each interval such that

    s_j(x) = a_j + b_j * (x-x_j) + c_j * (x-x_j)^2 + d_j * (x-x_j)^3,
    
    with a's being the first row, b's the second, and so on.
"""
function compute_natural_spline(xi, fi)
    np1 = length(xi)
    result = zeros(4, np1-1)
    # represent the spline on each interval as
    # ai + bi (x-xi) + ci (x-xi)^2 + di (x-xi)^2
    # result will hold 1 column per cell, with [ai; bi; ci; di]
    
    h = xi[2:np1]-xi[1:np1-1]
    
    # Let's build the linear system we derived as a tridiagonal system
    # Julia has a tool that just lets us build the diagonals.
    vd = zeros(np1)
    vl = zeros(np1-1)
    vu = zeros(np1-1)
    vd[1] = 1.0
    vd[end] = 1.0
    for j=2:np1-1
        vl[j] = h[j-1]
        vu[j] = h[j]
        vd[j] = 2*(h[j-1] + h[j])
    end
    M = LinearAlgebra.Tridiagonal(vl, vd, vu)

    # now we need to build the RHS
    rhs = zeros(np1)
    for j=2:np1-1
        rhs[j] = 3/h[j]*(fi[j+1]-fi[j]) - 3/h[j-1]*(fi[j]-fi[j-1])
    end
    
    # Solve for c values.
    cs = M \ rhs

    result[1, :] = fi[1:np1-1]  # a
    result[3, :] = cs[1:np1-1]  # c
    
    # b and d values
    for j=1:np1-2
        result[2,j] = (fi[j+1]-fi[j])/h[j] - h[j]*(2.0*result[3, j] + result[3, j+1])/3.0
        result[4,j] = (result[3,j+1]-result[3,j])/3.0/h[j]
    end
    j = np1-1
    result[2,j] = (fi[j+1]-fi[j])/h[j] - h[j]*(2.0*result[3, j])/3.0
    result[4,j] = -result[3,j]/3.0/h[j] 
    
    return result
end


"""
    evaluate_spline(xi, coeffs, x)

    Given the abscissae xi and the computed coeffs defining a
    cubic spline, evaluates the spline at a given point x.

    It first has to locate which interval x belongs to, and then
    evaluate the cubic polynomial on that interval.
"""
function evaluate_spline(xi, coeffs, x)
    n = length(xi)-1
    function find_i(xi, x)
        # This does a linear search.  Binary is better!
        for i=1:n
            if x >= xi[i] && x < xi[i+1]
                return i
            end
        end
        return n
    end
    i = find_i(xi, x)
    xx = x-xi[i]
    # This is Horner-factored for optimality!
    return xx*(xx*(xx*(coeffs[4,i]) + coeffs[3, i]) + coeffs[2, i]) + coeffs[1,i]
end


"""
    evaluate_spline_derivative(xi, coeffs, x)

    Given the abscissae xi and the computed coeffs defining a
    cubic spline, evaluates the derivative of the spline at a given point x.
    This could also be done by automatic differentiation on 
    the spline evaluation function.

    It first has to locate which interval x belongs to, and then
    evaluate the cubic polynomial on that interval.
"""
function evaluate_spline_derivative(xi, coeffs, x)
    n = length(xi)-1
    function find_i(xi, x)
        # This does a linear search.  Binary is better!
        for i=1:n
            if x >= xi[i] && x < xi[i+1]
                return i
            end
        end
        return n
    end
    i = find_i(xi, x)
    xx = x-xi[i]
    # This is Horner-factored for optimality!
    return xx*(xx*3*coeffs[4,i]+2*coeffs[3,i]) + coeffs[2,i]
end;


# xs = [0.1, 0.2, 0.3, 0.4];
# f(x) = x * cos(x) - 2 * x^2 + 3*x - 1;

# ys = f.(xs)

# spl = compute_natural_spline(xs, ys)

# using Plots

# plot(f, 0.1, 0.4)
# plot!(x->evaluate_spline(xs, spl, x), 0.1, 0.4)
# savefig("foo.png")
