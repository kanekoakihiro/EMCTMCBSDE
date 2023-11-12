abstract type AbstractGrid end



# Plane grid
struct Grid1D <: AbstractGrid
    grid
    N
    Nₗ
    Nᵣ
end



# Regular grid
abstract type AbstractInterval end
struct UnboundedInterval <: AbstractInterval end
struct LeftUnboundedInterval <: AbstractInterval end
struct RightUnboundedInterval <: AbstractInterval end
struct BoundedInterval <: AbstractInterval end

function AbstractInterval(domain::Vector)
    @assert length(domain) == 2
    @assert domain[1] < domain[2]
    if domain[1] == -Inf
        domain[2] == Inf ? (return UnboundedInterval()) : (return LeftUnboundedInterval())
    else
        domain[2] == Inf ? (return RightUnboundedInterval()) : (return BoundedInterval())
    end
end

function Nₗᵣ(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64})
    @assert center > lim[1] && center < lim[2]
    @assert Δₓ > 0
    @assert Nₓ > 0
    return Nₗᵣ(Δₓ, Nₓ, center, lim, AbstractInterval(lim))
end

function Nₗᵣ(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64}, ::UnboundedInterval)
    return (Nₗ = Nₓ, Nᵣ = Nₓ)
end

function Nₗᵣ(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64}, ::LeftUnboundedInterval)
    center+Nₓ*Δₓ > lim[2] ? Nᵣ = Int((lim[2]-center)÷Δₓ) : Nᵣ = Nₓ
    Nₗ = 2*Nₓ-Nᵣ
    return (Nₗ = Nₗ, Nᵣ = Nᵣ)
end

function Nₗᵣ(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64}, ::RightUnboundedInterval)
    center-Nₓ*Δₓ < lim[1] ? Nₗ = Int((center-lim[1])÷Δₓ) : Nₗ = Nₓ
    Nᵣ = 2*Nₓ-Nₗ
    return (Nₗ = Nₗ, Nᵣ = Nᵣ)
end

function Nₗᵣ(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64}, ::BoundedInterval)
    max_Nₗ = Int((center-lim[1])÷Δₓ)
    max_Nᵣ = Int((lim[2]-center)÷Δₓ)
    @assert max_Nₗ+max_Nᵣ >= 2*Nₓ "Nₓ or Δₓ need to be reduced."
    if max_Nₗ > max_Nᵣ
        Nᵣ = min(max_Nᵣ,Nₓ)
        Nₗ = 2*Nₓ - Nᵣ
    else
        Nₗ = min(max_Nₗ,Nₓ)
        Nᵣ = 2*Nₓ - Nₗ
    end
    return (Nₗ = Nₗ, Nᵣ = Nᵣ)
end

struct RegularGrid1D <: AbstractGrid
    grid
    N
    Nₗ
    Nᵣ
    function RegularGrid1D(Δₓ::Float64, Nₓ::Int64, center::Float64, lim::Vector{Float64})
        Nₗ, Nᵣ = Nₗᵣ(Δₓ, Nₓ, center, lim)
        grid = collect(range(center-Nₗ*Δₓ,center+Nᵣ*Δₓ,step=Δₓ))
        N = 2*Nₓ+1
        new(grid, N, Nₗ, Nᵣ)
    end
end



# Tavella-Randall grid
struct TavellaRandallGrid <: AbstractGrid
    grid
    N
    Nₗ
    Nᵣ
    function TavellaRandallGrid(g₁::Float64, g₂::Float64, left::Float64, center::Float64, right::Float64, Nₗ::Int64, Nᵣ::Int64)
        c₁ = asinh((left-center)/g₁)
        c₂ = asinh((center-right)/g₂)
        N = Nₗ+Nᵣ+1
        grid = zeros(N)
        for i=1:(Nₗ+1)
            grid[i] = center + g₁*sinh(c₁*(1-(i-1)/Nₗ))
        end
        for i=1:Nᵣ
            grid[i+Nₗ+1] = center - g₂*sinh(c₂*(i/Nᵣ))
        end
        new(grid, N, Nₗ, Nᵣ)
    end
end