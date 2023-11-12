using Interpolations

function sparse_indices(q::Int64, d::Int64)
    cartesian_indices = collect(Base.product([collect(1:(q+1)) for i in 1:d]...))
    res = filter(x->((sum(x)>=q-d+1)&&(sum(x)<=q)),cartesian_indices)
    if d == 1
        res = sparse_indices[1][1]
    end
    return res
end

function m(i)
    return (2^i)+1
end

function eq_nodes(i)
    return collect(2.0.*(0:(m(i)-1))./(m(i)-1.0).-1.0)
end

function rescale(x, scaled_center, left, right)
    y = x.+scaled_center .- sign.(x).*scaled_center.*x
    return (right-left)/2.0.*y.+(right+left)/2.0
end

struct SparseGridInterpolation
    Interpolant
    function SparseGridInterpolation(datalist, grids, spind, q, d, Nₜ, lefts, rights, centers)
        SG_coefs = (l->((-1)^(q-l))*binomial(d-1,q-l)).(sum.(spind))
        itps = []
        for t in 1:Nₜ+1
            itps_t = []
            for (grid, data) in zip(grids, datalist)
                num_points = [gr.N for gr in grid]
                evals = permutedims(reshape(data[:,t], tuple(reverse(num_points)...)), d:-1:1)
                itp = interpolate(tuple([gr.grid for gr in grid]...), evals, Gridded(Linear()))
                push!(itps_t, itp)
            end 
            push!(itps, itps_t)
        end
        Interpolant = (t, x)-> sum(SG_coefs.*[itp(x...) for itp in itps[t]])
        new(Interpolant)
    end
end