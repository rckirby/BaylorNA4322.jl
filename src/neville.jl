"""
    neville(xs, fs, x)

    Evaluates the Lagrange interpolant using nodes xs, values fs at a point x.
    Explicitly forms the entire Neville tableaux.
"""
function neville_table(xs, fs, x)
    @assert length(xs) == length(fs)

    r = length(xs)

    Q = zeros(r, r)
    Q[:, 1] = fs[:]
    for j=2:r
        for i=j:r
            Q[i, j] = ((x-xs[i-j+1]) * Q[i, j-1] - (x-xs[i]) * Q[i-1, j-1])/(xs[i] - xs[i-j+1])
        end
    end

    return Q[end, end]
end


"""
    neville(xs, fs, x)

    Evaluates the Lagrange interpolant using nodes xs, values fs at a point x.
"""
function neville(xs, fs, x)
    @assert length(xs) == length(fs)
    r = length(xs)
    Q = zeros(r)
    Q[:] = fs[:]

    # Note: this avoids the O(N^2) storage by overwriting previous computations
    # However, because of the data dependencies in the loop, we have to
    # go from bottom to top.
    # We can circumvent the O(N^2) storage but not work.
    for j=2:r
        for i=r:-1:j
            Q[i] = ((x-xs[i-j+1]) * Q[i] - (x-xs[i]) * Q[i-1])/(xs[i] - xs[i-j+1])
        end
    end

    return Q[end]
end
