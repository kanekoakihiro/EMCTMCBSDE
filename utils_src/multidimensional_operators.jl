using Kronecker
function grid_Kron(grids::Array{AbstractGrid, 1})
    d = length(grids)
    gr_ext = Array{Array{Float64,1},1}(undef, d)
    for dim in 1:d
        if dim == 1
            res = grids[1].grid ⊗ ones(prod([gr.N for gr in grids[2:d]]))
        else
            res = ones(prod([gr.N for gr in grids[1:dim-1]])) ⊗ grids[dim].grid
            if dim < d
                res = res ⊗ ones(prod([gr.N for gr in grids[dim+1:d]]))
            end
        end
        gr_ext[dim] = collect(res)[:,1]
    end
    return hcat(gr_ext...)
end

function D₁_Kron(grids::Array{AbstractGrid, 1}, D₁::Array{Array{Float64,2},1})
    d = length(grids)
    D₁_ext = Array{Array{Float64,2},1}(undef, d)
    for dim in 1:d
        if dim == 1
            res = D₁[1] ⊗ Matrix(I, prod([gr.N for gr in grids[2:d]]), prod([gr.N for gr in grids[2:d]]))
        else
            res = Matrix(I, prod([gr.N for gr in grids[1:dim-1]]), prod([gr.N for gr in grids[1:dim-1]])) ⊗ D₁[dim]
            if dim < d
                res = res ⊗ Matrix(I, prod([gr.N for gr in grids[dim+1:d]]), prod([gr.N for gr in grids[dim+1:d]]))
            end
        end
        D₁_ext[dim] = collect(res)
    end
    return D₁_ext
end

function D₂_Kron(grids::Array{AbstractGrid, 1}, D₁::Array{Array{Float64,2},1}, D₂::Array{Array{Float64,2},1})
    d = length(grids)
    D₂_ext = Array{Array{Array{Float64,2},1},1}(undef, d)
    for dim₁ in 1:d
        tmp = Array{Array{Float64,2},1}(undef, d)
        for dim₂ in 1:d
            if dim₁ == dim₂
                if dim₁ == 1
                    res = D₂[1] ⊗ Matrix(I, prod([gr.N for gr in grids[2:d]]), prod([gr.N for gr in grids[2:d]]))
                else
                    res = Matrix(I, prod([gr.N for gr in grids[1:dim₁-1]]), prod([gr.N for gr in grids[1:dim₁-1]])) ⊗ D₂[dim₁]
                    if dim₁ < d
                        res = res ⊗ Matrix(I, prod([gr.N for gr in grids[dim₁+1:d]]), prod([gr.N for gr in grids[dim₁+1:d]]))
                    end
                end
            else
                ldim = min(dim₁, dim₂)
                rdim = max(dim₁, dim₂)
                if ldim == 1
                    res = D₁[1]
                else
                    res = Matrix(I, prod([gr.N for gr in grids[1:ldim-1]]), prod([gr.N for gr in grids[1:ldim-1]])) ⊗  D₁[ldim]
                end
                if rdim-ldim == 1
                    res = res ⊗ D₁[rdim]
                else
                    res = res ⊗ Matrix(I, prod([gr.N for gr in grids[ldim+1:rdim-1]]), prod([gr.N for gr in grids[ldim+1:rdim-1]])) ⊗ D₁[rdim]
                end
                if rdim < d
                    res = res ⊗ Matrix(I, prod([gr.N for gr in grids[rdim+1:d]]), prod([gr.N for gr in grids[rdim+1:d]]))
                end
            end
            tmp[dim₂] = collect(res)
        end
        D₂_ext[dim₁] = tmp
    end
    return D₂_ext
end


function generateDiffAttributes(grids::Array{AbstractGrid, 1}, b::Function, σ::Function, f::Function, g::Function)
    d = length(grids)
    
    D₁_list = [(if gr.N > 1 collect(compute_1st_diff_operator(gr.grid)') else [0.0] end) for gr in grids]
    D₂_list = [(if gr.N > 1 collect(compute_2nd_diff_operator(gr.grid)') else [0.0] end) for gr in grids]    
    
    grid = grid_Kron(grids)
    N = size(grid)[1]
    𝐃₁ = D₁_Kron(grids, D₁_list)
    𝐃₂ = D₂_Kron(grids, D₁_list, D₂_list)
    
    ### Construction of the spatial grid for our CTMC.    
    𝐛 = b.(eachslice(grid, dims=1))
    𝛔 = σ.(eachslice(grid, dims=1))
    𝐚 = (mat->mat*mat').(𝛔)

    𝐛 = permutedims(hcat(𝐛...));
    𝛔 = permutedims(reshape(hcat(𝛔...), d, d, :), [3, 1, 2]);
    𝐚 = permutedims(reshape(hcat(𝐚...), d, d, :), [3, 1, 2]);

    𝐐 = zeros(N, N)
    for p in 1:d
        𝐐 += Diagonal(𝐛[:,p])*𝐃₁[p]
        for q in 1:d
            𝐐 += 0.5.*Diagonal(𝐚[:,p,q])*𝐃₂[p][q]
        end
    end
    
    𝛔ᵀ𝐃₁ = Array{Matrix{Float64}, 1}(undef, d)
    for p in 1:d
        mat = zeros(N, N)
        for q in 1:d
            mat += Diagonal(𝛔[:,q,p])*𝐃₁[q]
        end
        𝛔ᵀ𝐃₁[p] = mat
    end
    states = eachslice(grid, dims=1)
    Nonlinear(v, p, t) = f.(t, states, v, eachslice(hcat([mat*v for mat in 𝛔ᵀ𝐃₁]...),dims=1));
    Terminal = g.(states)
    return (N=prod([gr.N for gr in grids]), grid=grid, D₁=𝐃₁, D₂=𝐃₂, b=𝐛, σ=𝛔, a=𝐚, Q=𝐐, σD₁=𝛔ᵀ𝐃₁, Nonlinear=Nonlinear, Terminal=Terminal)
end