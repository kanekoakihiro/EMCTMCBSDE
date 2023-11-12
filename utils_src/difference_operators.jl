function compute_1st_diff_operator(grid::Vector{Float64})
    if length(grid) == 1
        return [0.0]
    end
    ΔG=diff(grid)
    @views coef1=@. -ΔG[2:end]/(ΔG[1:end-1]*(ΔG[1:end-1]+ΔG[2:end]))
    @views coef2=@. ΔG[1:end-1]/(ΔG[2:end]*(ΔG[1:end-1]+ΔG[2:end]))
    @views coef0=@. (ΔG[2:end]-ΔG[1:end-1])/(ΔG[1:end-1]*ΔG[2:end])
    
    dlower = cat([0], coef2, dims=1)
    d = cat([0],coef0,[0],dims=1)
    dupper = cat(coef1,[0],dims=1)
    return Tridiagonal(dlower, d, dupper)
end
compute_1st_diff_operator(grid::AbstractGrid) = compute_1st_diff_operator(grid.grid)

function compute_2nd_diff_operator(grid::Vector{Float64})
    if length(grid) == 1
        return [0.0]
    end
    ΔG=diff(grid)
    @views coef1=@. 2.0/(ΔG[1:end-1]*(ΔG[1:end-1]+ΔG[2:end]))
    @views coef2=@. 2.0/(ΔG[2:end]*(ΔG[1:end-1]+ΔG[2:end]))
    @views coef0=@. -2.0/(ΔG[1:end-1]*ΔG[2:end])
    
    dlower = cat([0], coef2, dims=1)
    d = cat([0],coef0,[0],dims=1)
    dupper = cat(coef1,[0],dims=1)
    return Tridiagonal(dlower, d, dupper)
end
compute_2nd_diff_operator(grid::AbstractGrid) = compute_2nd_diff_operator(grid.grid)

function generateDiffAttributes(grid::AbstractGrid, b::Function, σ::Function, f::Function, g::Function)
    𝐃₁ = compute_1st_diff_operator(grid.grid)'
    𝐃₂ = compute_2nd_diff_operator(grid.grid)'
    
    ### Construction of the spatial grid for our CTMC.
    𝐛 = b.(grid.grid)
    𝛔 = σ.(grid.grid)
    𝐚 = 𝛔.^2

    𝐐 = 𝐛.*𝐃₁ .+ (0.5).*𝐚.*𝐃₂
    𝛔𝐃₁ = 𝛔.*𝐃₁
    
    Nonlinear(v, p, t) = f.(t, grid.grid, v, 𝛔𝐃₁*v)
    Terminal = g.(grid.grid)
    return (N=grid.N, grid=grid.grid, D₁=𝐃₁, D₂=𝐃₂, b=𝐛, σ=𝛔, a=𝐚, Q=𝐐, σD₁=𝛔𝐃₁, Nonlinear=Nonlinear, Terminal=Terminal)
end