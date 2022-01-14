function neville_table(xs, fs, x)
    @assert length(xs) == length(fs)

    # one more than the polynomial degree!
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

# Note: this avoids the O(N^2) storage by overwriting previous computations
# However, because of the data dependencies in the loop, we have to
# go from bottom to top.
function neville(xs, fs, x)
    @assert length(xs) == length(fs)
    r = length(xs)
    Q = zeros(r)
    Q[:] = fs[:]
    for j=2:r
        for i=r:-1:j
            Q[i] = ((x-xs[i-j+1]) * Q[i] - (x-xs[i]) * Q[i-1])/(xs[i] - xs[i-j+1])
        end
    end

    return Q[end]
end
    

    
